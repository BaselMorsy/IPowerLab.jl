include("meta_solver.jl")

function build_DOPF_MP!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    return build_full_DOPF_model!(grid,SimulationSettings,prerequisites_data)
end

function build_DOPF_SP!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    solved_MP_model::Model, t::Int, k::Int)
    prerequisites_data_instance = deepcopy(prerequisites_data)
    prerequisites_data_instance.time_horizon = [t]
    prerequisites_data_instance.k_t[t] = [k]

    # Fixing schedules and emptying commitable/non-commitable gen lists:
    # p_gen_ac = JuMP.value.(solved_MP_model[:p_gen_ac][:,1,t])
    if k ∉ prerequisites_data_instance.contingency_redispatch
        fixed_schedules = Dict()
        for g in prerequisites_data_instance.ac_gen_ids
            push!(fixed_schedules, g => Dict(t => JuMP.value(solved_MP_model[:p_gen_ac][g,1,t])))
        end
    end
    prerequisites_data_instance.fixed_schedules = fixed_schedules
    prerequisites_data_instance.commitable_gen_ids = []
    prerequisites_data_instance.non_commitable_gen_ids = []
    prerequisites_data_instance.fixed_commitments = Dict()
    
    # Fixing relevant topology:
    fixed_topology = Dict()

    if SimulationSettings.transmission_switching == [:pre] && length(prerequisites_data_instance.ac_active_dynamic_branch_ids) != 0
        prerequisites_data_instance.ac_fixed_dynamic_branch_ids = prerequisites_data_instance.ac_active_dynamic_branch_ids
        prerequisites_data_instance.ac_active_dynamic_branch_ids = []
        for branch_id in prerequisites_data_instance.ac_fixed_dynamic_branch_ids
            push!(fixed_topology, branch_id => Dict(k => Dict(t => JuMP.value(solved_MP_model[:z_l][branch_id, 1, t]))))
        end
    end

    if SimulationSettings.substation_switching["splitting"] == [:pre] && length(prerequisites_data_instance.ac_active_coupler_ids) != 0
        prerequisites_data_instance.ac_fixed_coupler_ids = prerequisites_data_instance.ac_active_coupler_ids
        prerequisites_data_instance.ac_active_coupler_ids = []
        for branch_id in prerequisites_data_instance.ac_fixed_coupler_ids
            push!(fixed_topology, branch_id => Dict(k => Dict(t => JuMP.value(solved_MP_model[:z_c][branch_id, 1, t]))))
        end
    end

    if SimulationSettings.substation_switching["reconf"] == [:pre] && length(prerequisites_data_instance.ac_active_reconf_ids) != 0
        prerequisites_data_instance.ac_fixed_reconf_ids = prerequisites_data_instance.ac_active_reconf_ids
        prerequisites_data_instance.ac_active_reconf_ids = []
        for branch_id in prerequisites_data_instance.ac_fixed_reconf_ids
            push!(fixed_topology, branch_id => Dict(k => Dict(t => JuMP.value(solved_MP_model[:z_r][branch_id, 1, t]))))
        end
    end

    prerequisites_data_instance.fixed_topology = fixed_topology

    # Building model:
    return build_snapshot_DOPF_model!(solved_MP_model, grid, SimulationSettings, prerequisites_data_instance, k, t)
end

function topology_warm_start!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    solved_MP_model::Model)

    starting_topology = Dict()
    S0 = Dict()
    S1 = Dict()

    if SimulationSettings.substation_switching["reconf"] == [:pre]
        for t in prerequisites_data.time_horizon
            for branch_id in prerequisites_data.ac_active_reconf_ids
                push!(starting_topology, branch_id => Dict(1 => Dict(t => JuMP.value(solved_MP_model[:z_r][branch_id, 1, t]))))
                if JuMP.value(solved_MP_model[:z_r][branch_id, 1, t]) == 0
                    if t ∉ keys(S0)
                        push!(S0, t => Dict())
                        if 1 ∉ keys(S0[t])
                            push!(S0[t], 1 => [branch_id])
                        end
                    else
                        if 1 ∉ keys(S0[t])
                            push!(S0[t], 1 => [branch_id])
                        else
                            push!(S0[t][1], branch_id)
                        end
                    end
                else
                    if t ∉ keys(S1)
                        push!(S1, t => Dict())
                        if 1 ∉ keys(S1[t])
                            push!(S1[t], 1 => [branch_id])
                        end
                    else
                        if 1 ∉ keys(S1[t])
                            push!(S1[t], 1 => [branch_id])
                        else
                            push!(S1[t][1], branch_id)
                        end
                    end
                end
            end
        end
    elseif SimulationSettings.substation_switching["reconf"] == [:pre,:post]
        for t in prerequisites_data.time_horizon
            for k in prerequisites_data.k_t[t]
                for branch_id in prerequisites_data.ac_active_reconf_ids
                    push!(starting_topology, branch_id => Dict(k => Dict(t => JuMP.value(solved_MP_model[:z_r][branch_id, 1, t]))))
                    if JuMP.value(solved_MP_model[:z_r][branch_id, k, t]) == 0
                        if t ∉ keys(S0)
                            push!(S0, t => Dict())
                            if k ∉ keys(S0[t])
                                push!(S0[t], k => [branch_id])
                            end
                        else
                            if k ∉ keys(S0[t])
                                push!(S0[t], k => [branch_id])
                            else
                                push!(S0[t][k], branch_id)
                            end
                        end
                    else
                        if t ∉ keys(S1)
                            push!(S1, t => Dict())
                            if k ∉ keys(S1[t])
                                push!(S1[t], k => [branch_id])
                            end
                        else
                            if k ∉ keys(S1[t])
                                push!(S1[t], k => [branch_id])
                            else
                                push!(S1[t][k], branch_id)
                            end
                        end
                    end
                end
            end
        end
    end

    return starting_topology, S0, S1
end

function solve_decomposed_model!(model::Model; t_max=Inf)
    JuMP.set_attribute(model, "TimeLimit", t_max)
    optimize!(model)
    return model
end

function fix_converter_flows!( prerequisites_data::DOPF_Prerequisites, solved_MP_model::Model, SP_model::Model, t::Int, k::Int)
   JuMP.@constraint(SP_model, fix_converter_flows[g in prerequisites_data.conv_ac_side_virtual_gen_ids],
        SP_model[:p_conv_ac][g,k,t] == JuMP.value(solved_MP_model[:p_conv_ac][g,1,t])) 
end

function solve_DOPF_CCG!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false, ϵ=0.1, update_order_book=false)

    δ = Inf
    
    iter_count = 0
    T = SimulationSettings.time_horizon
    K_all = deepcopy(prerequisites_data.k)

    t_master = []
    t_sub = []

    δ_it = []
    LS_it = [] # to be calculated later

    solved_MP_model = []

    MP_model = []
    MP_models = Dict()
    SP_models = Dict()
    k_last = Dict()
    starting_topology = Dict()
    S0 = []
    S1 = []
    S0_all = []
    S1_all = []
    τ_visited = []
    τ_cut = []
    reconf_lines_set_to_1 = []

    LB_it = []
    UB_it = []

    if get(SimulationSettings.Meta_solver_instance.misc,"k_final",[]) != []
        prerequisites_data.k_t = SimulationSettings.Meta_solver_instance.misc["k_final"]
    end

    if get(SimulationSettings.Meta_solver_instance.misc,"starting_topology",[]) != []
        starting_topology = SimulationSettings.Meta_solver_instance.misc["starting_topology"]
    end

    if SimulationSettings.converter_modularization == :continuous && length(grid.N_conv_duplets) > 0
        Converter_Duplets = grid.Converter_Duplets
        for pair_id in keys(Converter_Duplets)
            main_ac_bus = grid.Converters[Converter_Duplets[pair_id][1]].AC_Bus_ID
            if main_ac_bus ∈ keys(grid.Substations)
                module_1 = grid.Converters[Converter_Duplets[pair_id][1]]
                module_2 = grid.Converters[Converter_Duplets[pair_id][2]]
                g_bus_1 = grid.Generators[module_1.gen_ac_id].GenBus_ID
                g_bus_2 = grid.Generators[module_2.gen_ac_id].GenBus_ID
                line_1 = grid.Buses[g_bus_1].ConnectedLinesIDs[1]
                line_2 = grid.Buses[g_bus_2].ConnectedLinesIDs[2]
                @assert grid.Branches[line_1].Fr_bus_ID != grid.Branches[line_2].Fr_bus_ID

                push!(reconf_lines_set_to_1, line_1)
                push!(reconf_lines_set_to_1, line_2)
            end
        end
    end
    #Main loop
    while δ ≥ ϵ
        iter_count += 1

        println("Debug msg @ iteration " * string(iter_count))
        # Solve master-problem
        MP_model = build_DOPF_MP!(grid, SimulationSettings, prerequisites_data)
        
        ## apply starting topology
        if length(keys(starting_topology)) > 0
            for t in prerequisites_data.time_horizon
                for r in prerequisites_data.ac_active_reconf_ids
                    JuMP.set_start_value(MP_model[:z_r][r,1,t], starting_topology[r][1][t])
                end
                # for c in prerequisites_data.ac_active_coupler_ids
                #     for k in prerequisites_data.k_t[t][2:end]
                #         if k in k_last[(iter_count-1,t)]
                #             value = JuMP.value(MP_models[iter_count-1][:z_c][c,k,t])
                #             JuMP.set_start_value(MP_model[:z_c][c,k,t], value)
                #         # else
                #         #     value = JuMP.value(SP_models[(iter_count-1,t,k)][:z_c][c,k,t])
                #         #     JuMP.set_start_value(MP_model[:z_c][c,k,t], value)
                #         end
                        
                #     end
                # end
            end
            # if iter_count > 2
            #     if SimulationSettings.substation_switching["reconf"] == [:pre]
            #         JuMP.@constraint(MP_model, explored_topologies[i in 1:length(keys(S0_all))],
            #             sum([MP_model[:z_r][r,1,t] for t in prerequisites_data.time_horizon, r in prerequisites_data.ac_active_reconf_ids if r ∈ S0_all[i][t][1]])
            #             + sum([1-MP_model[:z_r][r,1,t] for t in prerequisites_data.time_horizon, r in prerequisites_data.ac_active_reconf_ids if r ∈ S1_all[i][t][1]]) ≥ 1)
            #     elseif SimulationSettings.substation_switching["reconf"] == [:pre,:post]
            #         JuMP.@constraint(MP_model, explored_topologies[i in 1:length(keys(S0_all))],
            #             sum([MP_model[:z_r][r,k,t] for t in prerequisites_data.time_horizon, k in prerequisites_data.k if k ∈ prerequisites_data.k_t[t] && r ∈ S0_all[i][t][k]])
            #             + sum([1-MP_model[:z_r][r,k,t] for t in prerequisites_data.time_horizon, k in prerequisites_data.k if k ∈ prerequisites_data.k_t[t] && r ∈ S1_all[i][t][k]]) ≥ 1)
            #     end
            # end
        end
        # if iter_count > 1
        #     # apply lower bound on cost
        #     GenBids = prerequisites_data.Order_Book.Gen_bids
        #     LoadBids = prerequisites_data.Order_Book.Load_bids
        #     JuMP.@constraint(MP_model, cost_lower_bound, sum([MP_model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in prerequisites_data.non_commitable_gen_ids, t in prerequisites_data.time_horizon], init=0)
        #         + sum([MP_model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] + grid.Generators[g].C0*MP_model[:u_gt][g,t]+MP_model[:α_gt][g,t]*grid.Generators[g].start_up_cost+MP_model[:β_gt][g,t]*grid.Generators[g].shut_down_cost for g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], init=0)
        #         + sum([MP_model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon], init=0)
        #         + sum([MP_model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in keys(prerequisites_data.fixed_schedules), t in prerequisites_data.time_horizon], init=0)
        #         + sum([MP_model[:p_gen_dc][g,1,t]*grid.DCGenerators[g].C1 for g in prerequisites_data.dc_gen_ids, t in prerequisites_data.time_horizon], init=0)
        #         + sum([MP_model[:Γ][t] for t in prerequisites_data.time_horizon], init=0) ≥ sum(LB_it[end], init=0))
        # end

        if reconf_lines_set_to_1 != []
            for t in prerequisites_data.time_horizon
                for r in reconf_lines_set_to_1
                    JuMP.fix(MP_model[:z_r][r,1,t], 1; force = true)
                end
            end
        end
        t_master_now = @elapsed solved_MP_model = solve_decomposed_model!(MP_model; t_max = get(SimulationSettings.Meta_solver_instance.misc,["t_max"],Inf))
        
        if !JuMP.has_values(solved_MP_model)
            status = process_last_MP!(grid, solved_MP_model, prerequisites_data, SimulationSettings, order_book, update_grid, update_order_book)
            return solved_MP_model, status
        end
        push!(t_master, t_master_now)
        push!(MP_models, iter_count => solved_MP_model)

        LB_t = [calculate_opex_t(solved_MP_model,grid, prerequisites_data, t; include_load_shedding=true, include_commitment_cost=false) for t in T]

        push!(LB_it, LB_t)

        UB_tk = zeros(length(T), length(K_all))

        if SimulationSettings.Parallels
            t_sub_now = @elapsed Threads.@threads for t in T
                Threads.@threads for k in sort(collect(setdiff(Set(K_all),Set(prerequisites_data.k_t[t])))) 
                    SP_model = build_DOPF_SP!(grid, SimulationSettings, prerequisites_data, solved_MP_model, t, k)
                    if !SimulationSettings.dynamic_converter_control
                        fix_converter_flows!(prerequisites_data, solved_MP_model, SP_model, t, k)
                    end
                    solved_SP_model = solve_decomposed_model!(SP_model)
                    push!(SP_models, (iter_count,t,k) => solved_SP_model)
                    UB_tk[t,k] = JuMP.objective_value(solved_SP_model)
                end
            end
        else
            t_sub_now = @elapsed for t in T
                for k in sort(collect(setdiff(Set(K_all),Set(prerequisites_data.k_t[t])))) 
                    SP_model = build_DOPF_SP!(grid, SimulationSettings, prerequisites_data, solved_MP_model, t, k)
                    if !SimulationSettings.dynamic_converter_control
                        fix_converter_flows!(prerequisites_data, solved_MP_model, SP_model, t, k)
                    end
                    solved_SP_model = solve_decomposed_model!(SP_model)
                    push!(SP_models, (iter_count,t,k) => solved_SP_model)
                    UB_tk[t,k] = JuMP.objective_value(solved_SP_model)
                end
            end
        end
        push!(t_sub, t_sub_now)
        UB_t = zeros(length(T),1)
        
        for t in T
            # 1 should be a data structure ::CCG_Settings which can include nk[t] (number of k's to include) or a DRL agent that chooses the k's
            UB_t[t] = maximum(UB_tk[t,:])
        end

        push!(UB_it, UB_t)

        δ = maximum((UB_t - LB_t) ./ LB_t)

        if δ ≥ ϵ
            for t in T
                next_k = find_next_k(UB_tk[t,:],prerequisites_data.k_t[t][2:end], 1)
                println("Next k for t = "*string(t)*": "*string(next_k))
                println("Current k's for t = "*string(t)*": "*string(prerequisites_data.k_t[t]))
                push!(k_last, (iter_count,t) => prerequisites_data.k_t[t])
                append!(prerequisites_data.k_t[t], next_k)
            end
            println(LB_t)
            println(UB_t)
            starting_topology, S0, S1 = topology_warm_start!(grid, SimulationSettings, prerequisites_data, solved_MP_model)
            push!(τ_visited, starting_topology)
            if iter_count > 1
                if starting_topology == τ_visited[end-1] # meaning that the same topology is still useful
                    # Do nothing because we are not sure if this topology will fail or not in the next iteration
                else
                    # the topology now is different from the last iterations's topology -> meaning that τ_{i-1} is obsolete now
                    push!(S0_all, S0)
                    push!(S1_all, S1)
                    push!(τ_cut, τ_visited[end-1])
                end
            end
        end
        push!(δ_it, δ)
        
        println("==================================")
        println("|Iteration Count|δ|")
        println(iter_count, "  ", δ)
        println("==================================")
    end

    meta_solver = MetaSolver(MP_models=MP_models[iter_count],t_master=t_master, t_sub=t_sub, δ_it=δ_it,
        LB_it=LB_it, UB_it=UB_it, misc=Dict("k_final"=>prerequisites_data.k_t, "starting_topology"=>starting_topology))

    SimulationSettings.Meta_solver_instance = meta_solver 
    status = process_last_MP!(grid, solved_MP_model, prerequisites_data, SimulationSettings, order_book, update_grid, update_order_book)
    return solved_MP_model, status
end

function process_last_MP!(grid, model, prerequisites_data, SimulationSettings, order_book, update_grid, update_order_book)
    
    if JuMP.has_values(model)
        push!(grid.Operating_Cost, JuMP.objective_value(model))
        flag = haskey(model, :p_ls_ac)
        if update_order_book
            for t in prerequisites_data.time_horizon
                for g in prerequisites_data.Order_Book.Gen_ids
                    prerequisites_data.Order_Book.Schedule["gen"][g][t] = JuMP.value.(model[:p_gen_ac])[g,1,t]
                end

                if flag
                    for d in prerequisites_data.ac_load_shedding_ids
                        prerequisites_data.Order_Book.Schedule["load"][d][t] = -JuMP.value.(model[:p_ls_ac][d,1,t])
                    end
                else
                    for d in prerequisites_data.ac_load_shedding_ids
                        prerequisites_data.Order_Book.Schedule["load"][d][t] = 0
                    end
                end  
            end
        end

        if update_grid
            DOPF_post_processing!(model, grid, SimulationSettings, prerequisites_data, order_book)
        end
        return true
    else
        if update_order_book
            for t in prerequisites_data.time_horizon
                for g in prerequisites_data.Order_Book.Gen_ids
                    prerequisites_data.Order_Book.Schedule["gen"][g][t] = 0
                end


                for d in prerequisites_data.ac_load_shedding_ids
                    prerequisites_data.Order_Book.Schedule["load"][d][t] = 0
                end
            end
        end
        push!(grid.Operating_Cost, -1)
        return false
    end
end

function find_next_k(my_vector::Vector,my_other_vector::Vector, n::Int)
    # Create a filtered vector of indices for elements not in my_other_vector
    filtered_indices = filter(i -> !(my_vector[i] in my_other_vector), 1:length(my_vector))
    
    # Check if n is greater than the size of the filtered vector
    if n > length(filtered_indices)
        return filtered_indices  # Return all indices if n is larger
    end
    
    # Sort the filtered indices based on my_vector and return the largest n indices
    sorted_indices = sortperm(my_vector[filtered_indices], rev=true)
    largest_n_indices = filtered_indices[sorted_indices[1:n]]
    
    return largest_n_indices
end
