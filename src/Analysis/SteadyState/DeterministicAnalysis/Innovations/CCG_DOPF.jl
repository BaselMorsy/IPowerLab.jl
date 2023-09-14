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
    return build_snapshot_DOPF_model!(grid, SimulationSettings, prerequisites_data_instance, k=k)
end

function solve_decomposed_model!(model::Model)
    optimize!(model)
    return model
end

function solve_DOPF_CCG!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false, ϵ=0.1, update_order_book=false)

    δ = Inf
    
    iter_count = 0
    T = SimulationSettings.time_horizon
    K_all = prerequisites_data.k

    t_master = []
    t_sub = []

    δ_it = []
    LS_it = [] # to be calculated later

    solved_MP_model = []

    MP_model = []
    SP_models = Dict()
    #Main loop
    while δ ≥ ϵ
        iter_count = iter_count + 1

        println("Debug msg @ iteration " * string(iter_count))
        # Solve master-problem

        MP_model = build_DOPF_MP!(grid, SimulationSettings, prerequisites_data)
        t_master_now = @elapsed solved_MP_model = solve_decomposed_model!(MP_model)
        push!(t_master, t_master_now)

        LB_t = [calculate_opex_t(solved_MP_model,grid, prerequisites_data, t; include_load_shedding=false, include_commitment_cost=false) for t in T]

        UB_tk = zeros(length(T), length(K_all))

        if SimulationSettings.Parallels
            t_sub_now = @elapsed Threads.@threads for t in T
                Threads.@threads for k in collect(setdiff(Set(K_all),Set(prerequisites_data.k_t[t])))
                    SP_model = build_DOPF_SP!(grid, SimulationSettings, prerequisites_data, solved_MP_model, t, k)
                    solved_SP_model = solve_decomposed_model!(SP_model)
                    push!(SP_models, (t,k) => solved_SP_model)
                    UB_tk[t,k] = JuMP.objective_value(solved_SP_model)
                end
            end
        else
            t_sub_now = @elapsed for t in T
                for k in collect(setdiff(Set(K_all),Set(prerequisites_data.k_t[t])))
                    SP_model = build_DOPF_SP!(grid, SimulationSettings, prerequisites_data, solved_MP_model, t, k)
                    solved_SP_model = solve_decomposed_model!(SP_model)
                    push!(SP_models, (t,k) => solved_SP_model)
                    UB_tk[t,k] = JuMP.objective_value(solved_SP_model)
                end
            end
        end
        push!(t_sub, t_sub_now)
        UB_t = zeros(length(T),1)
        
        for t in T
            # 1 should be a data structure ::CCG_Settings which can include nk[t] (number of k's to include) or a DRL agent that chooses the k's
            append!(prerequisites_data.k_t[t], find_next_k(UB_tk[t,:], 1))
            UB_t[t] = maximum(UB_tk[t,:])
        end

        δ = maximum(abs.(UB_t - LB_t) ./ LB_t)
        push!(δ_it, δ)
        
        println("|Iteration Count|δ|")
        println(iter_count, "  ", δ)
    end
    status = process_last_MP!(model, update_grid, update_order_book)
    return solved_MP_model, status
end

function process_last_MP!(model, update_grid, update_order_book)
    
    if JuMP.has_values(model)
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
            push!(grid.Operating_Cost, JuMP.objective_value(model))
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

        if update_grid
            push!(grid.Operating_Cost, -1)
        end

        return false
    end
end

function find_next_k(v::Vector, nk::Int)
    return sortperm(v, rev=true)[1:nk]
end
