
@with_kw mutable struct Unit_Commitment_Simulation_Settings

    T = 24
    reconf = collect(1:24)
    splitting = collect(1:24)
    transmission_switching = []    
end

function solve_NCUC!(grid ::PowerGrid;T_=nothing,reference_gen=nothing)
    
    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]

    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    normally_opened_reconf_lines_per_coupler = Dict()

    for substation_id in keys(grid.Substations)
        substation_obj = grid.Substations[substation_id]
        coupler_id = substation_obj.Reconf_CouplerLines_IDs[1]
        reconf_ids = substation_obj.Reconf_AuxLines_IDs
        normally_opened_ids = []
        for reconf_id in reconf_ids
            if grid.Branches[reconf_id].GeneralSwitch.SwitchingStatus == 0
                push!(normally_opened_ids)
            end
        end
        push!(normally_opened_reconf_lines_per_coupler, coupler_id => normally_opened_ids)
    end

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", 900*3)
    # set_optimizer_attribute(m, "MIPGap", 0.01)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set, t in T] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, p[g in Gen_set,t in T])
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, t in T; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, u_gt[g in Gen_set, t in T], Bin)
    JuMP.@variable(m, α_gt[g in Gen_set, t in T], Bin) # start-up
    JuMP.@variable(m, β_gt[g in Gen_set, t in T], Bin) # shut-down
    JuMP.@variable(m, z_c[coupler = Coupler_set, t in T], Bin)
    JuMP.@variable(m, z_l[reconf = Reconf_set, t in T], Bin)

    JuMP.@constraint(m, Pg_max[g in Gen_set, t in T], p[g,t] ≤ grid.Generators[g].Pg_max*u_gt[g,t])
    JuMP.@constraint(m, Pg_min[g in Gen_set, t in T], p[g,t] ≥ grid.Generators[g].Pg_min*u_gt[g,t])
    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    JuMP.@constraint(m, logic[g in Gen_set, t in T[2:T_]], u_gt[g,t]-u_gt[g,t-1] == α_gt[g,t] - β_gt[g,t])
    JuMP.@constraint(m, logic_init[g in Gen_set], u_gt[g,1] == α_gt[g,1] + β_gt[g,1])
    
    JuMP.@constraint(m, MUT[g in Gen_set, t in grid.Generators[g].min_up_time:T_],
        sum(α_gt[g,i] for i in t-grid.Generators[g].min_up_time+1:t) ≤ u_gt[g,t])

    JuMP.@constraint(m, MDT[g in Gen_set, t in grid.Generators[g].min_down_time:T_],
        sum(β_gt[g,i] for i in t-grid.Generators[g].min_down_time+1:t) ≤ 1-u_gt[g,t])
    # 2. Constraints

    for t in T
        for coupler in Coupler_set
            set_start_value(z_c[coupler, t], grid.Branches[coupler].GeneralSwitch.SwitchingStatus_t[t])
        end
        for reconf in Reconf_set
            set_start_value(z_l[reconf, t ], grid.Branches[reconf].GeneralSwitch.SwitchingStatus_t[t])
        end
    end

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[t in T], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,t] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set, t in T],
        sum(pij[i,j,t] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        pij[i,j,t] == Sbase*(B[i,j])*(δ[i,t]-δ[j,t]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Branch_nodes], 
        pij[i,j,t] == -pij[j,i,t])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,t] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set, t in T],sum(z_l[l,t] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
    
    M = 1e3 # BigM number
    JuMP.@constraint(m,switched_off_reconf[coupler in Coupler_set, t in T], sum(z_l[l,t] for l in normally_opened_reconf_lines_per_coupler[coupler]) ≤ (1-z_c[coupler, t])*M)

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[j,t]-(1-z_l[Reconf_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_l[Reconf_dict[ Set([i,j])],t])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[j,t] - (1-z_c[Coupler_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_c[Coupler_dict[ Set([i,j])],t]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        -z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        pij[i,j,t] ≤ z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        pij[i,j,t] ≤ z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*u_gt[g,t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))

    JuMP.optimize!(m)

    solution_status = raw_status(m)

    if has_values(m)
        grid.Operating_Cost = JuMP.objective_value(m)
        Pgt = JuMP.value.(p)
        Pijt = JuMP.value.(pij)
        δt = JuMP.value.(δ)
        Z_c = JuMP.value.(z_c)
        Z_l = JuMP.value.(z_l)
        u_gt = JuMP.value.(u_gt)

        for g in keys(grid.Generators)
            grid.Generators[g].GeneralSwitch.SwitchingStatus_t = [u_gt[g,t] for t in T]
            grid.Generators[g].Pg_t = [Pgt[g,t] for t in T]
            grid.Generators[g].Qg_t = zeros(1,T_)
        end

        for b in keys(grid.Branches)
            Fr_bus = grid.Branches[b].Fr_bus_ID
            To_bus = grid.Branches[b].To_bus_ID
            grid.Branches[b].PowerFlow_ij_t = [Pijt[Fr_bus,To_bus,t] for t in T]
            grid.Branches[b].PowerFlow_ji_t = [Pijt[To_bus,Fr_bus,t] for t in T]
            grid.Branches[b].losses_P_t = [Pijt[Fr_bus,To_bus,t]+Pijt[To_bus,Fr_bus,t] for t in T]
            
            grid.Branches[b].ReactFlow_ij_t = zeros(1,T_)
            grid.Branches[b].ReactFlow_ji_t = zeros(1,T_)
            grid.Branches[b].losses_Q_t = zeros(1,T_)

            branch_type = grid.Branches[b].BranchType
            if branch_type == 1 # reconf line
                grid.Branches[b].GeneralSwitch.SwitchingStatus_t = [Z_l[b,t] for t in T]
            elseif branch_type == 2 # coupler line
                grid.Branches[b].GeneralSwitch.SwitchingStatus_t = [Z_c[b,t] for t in T]
            end
        end

        for bus in keys(grid.Buses)
            grid.Buses[bus].δ_t = [δt[bus,t] for t in T]
            grid.Buses[bus].V_magnitude_t = ones(1,T_)
        end
    end
    
    return solution_status
    
end

function solve_UC_gridless!(grid ::PowerGrid; T_=nothing)

    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_
    Gen_set = keys(grid.Generators)

    m = JuMP.Model(Gurobi.Optimizer)

    # 1. Variables
    JuMP.@variable(m, p[g in Gen_set,t in T])
    JuMP.@variable(m, u_gt[g in Gen_set, t in T], Bin)
    JuMP.@variable(m, α_gt[g in Gen_set, t in T], Bin) # start-up
    JuMP.@variable(m, β_gt[g in Gen_set, t in T], Bin) # shut-down

    JuMP.@constraint(m, Pg_max[g in Gen_set, t in T], p[g,t] ≤ grid.Generators[g].Pg_max*u_gt[g,t])
    JuMP.@constraint(m, Pg_min[g in Gen_set, t in T], p[g,t] ≥ grid.Generators[g].Pg_min*u_gt[g,t])
    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    JuMP.@constraint(m, logic[g in Gen_set, t in T[2:T_]], u_gt[g,t]-u_gt[g,t-1] == α_gt[g,t] - β_gt[g,t])
    JuMP.@constraint(m, logic_init[g in Gen_set], u_gt[g,1] == α_gt[g,1] + β_gt[g,1])
    
    JuMP.@constraint(m, MUT[g in Gen_set, t in grid.Generators[g].min_up_time:T_],
        sum(α_gt[g,i] for i in t-grid.Generators[g].min_up_time+1:t) ≤ u_gt[g,t])

    JuMP.@constraint(m, MDT[g in Gen_set, t in grid.Generators[g].min_down_time:T_],
        sum(β_gt[g,i] for i in t-grid.Generators[g].min_down_time+1:t) ≤ 1-u_gt[g,t])
    # 2. Constraints

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m,power_balance[t in T], sum(p[g,t] for g in Gen_set) == sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads)))
    
    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*u_gt[g,t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))

    JuMP.optimize!(m)

    solution_status = raw_status(m)

    if has_values(m)
        grid.Operating_Cost = JuMP.objective_value(m)
        Pgt = JuMP.value.(p)
        u_gt = JuMP.value.(u_gt)

        for g in keys(grid.Generators)
            grid.Generators[g].GeneralSwitch.SwitchingStatus_t = [u_gt[g,t] for t in T]
            grid.Generators[g].Pg_t = [Pgt[g,t] for t in T]
            grid.Generators[g].Qg_t = zeros(1,T_)
        end
    end
    
    return solution_status
end

function solve_UC_temporal_decomposition!(grid ::PowerGrid, SP_method;T_=nothing,reference_gen=nothing)
    
    T_congested_visited = Set([])
    t_reconf_set = Set([])
    ϵ = 1e-4
    UB = 9999999
    LB = -UB
    LB_arr = []
    UB_arr = []
    while abs((UB-LB)/LB) > ϵ
        
        UB, T_congested, t_l_congested_dict = solve_UC_MP!(grid,t_reconf_set,T_=T_,reference_gen=reference_gen)
        T_congested_new = setdiff(T_congested,T_congested_visited)
        t_l_congested_dict_reduced = Dict()
        for t in collect(T_congested_new)
            t_l_congested_dict_reduced[t] = t_l_congested_dict[t]
        end
        println(collect(T_congested_visited))
        println(collect(T_congested))
        println(collect(T_congested_new))
        push!(UB_arr,UB)
        if isempty(collect(T_congested_new))
            println("A7A")
            return UB
        end
        if SP_method == "OBS"
            LB, t_reconf_new = solve_UC_SP_OBS!(grid,T_congested_new,UB,reference_gen=reference_gen)
        elseif SP_method == "MC"
            LB, t_reconf_new = solve_UC_SP_MC!(grid, T_congested_new, t_l_congested_dict_reduced, UB, reference_gen=reference_gen)
        elseif SP_method == "ED"
            LB, t_reconf_new = solve_UC_SP_ED!(grid, T_congested_new, UB)
        else
            println("Invalid sub-problem method")
        end
        println(t_reconf_new)
        println((UB,LB))
        push!(LB_arr,LB)
        push!(t_reconf_set,t_reconf_new)
        push!(T_congested_visited,t_reconf_new)
    end

    return UB
end

function solve_NCUC_multi_level!(grid ::PowerGrid;ϵ=1e-6,time_limit=600,T_=nothing,reference_gen=nothing)
    UB = 999999
    LB = -999999
    
    LB_arr = []
    UB_arr = []
    while abs((UB-LB)/LB) > ϵ
        
        UB, α_gt, β_gt = solve_NCUC_L1!(grid,T_=T_,reference_gen=reference_gen)
        LB = solve_NCUC_L2_decomposed!(grid,α_gt, β_gt,time_limit,T_=T_,reference_gen=reference_gen)
        
        println((UB,LB))
        push!(UB_arr,UB)
        push!(LB_arr,LB)
    end

    return LB, UB_arr, LB_arr
end

function solve_NCUC_bilevel!(grid ::PowerGrid;ϵ=1e-6,time_limit=600,T_=nothing,reference_gen=nothing)
    UB = 999999
    LB = -999999
    
    LB_arr = []
    UB_arr = []
    while abs((UB-LB)/LB) > ϵ
        
        UB, α_gt, β_gt = solve_NCUC_L1!(grid,T_=T_,reference_gen=reference_gen)
        LB = solve_NCUC_L2!(grid,α_gt, β_gt,time_limit,T_=T_,reference_gen=reference_gen)
        
        println((UB,LB))
        push!(UB_arr,UB)
        push!(LB_arr,LB)
    end

    return LB, UB_arr, LB_arr
end

function solve_NCUC_L2_decomposed!(grid ::PowerGrid,α_gt,β_gt, time_limit;T_=nothing,reference_gen=nothing)
    substation_configs = Dict()
    for substation in keys(grid.Substations)
        z_l,z_c = solve_NCUC_L2_SP(grid,α_gt,β_gt,substation,time_limit, T_ = T_, reference_gen=reference_gen)
        config = Dict()
        push!(config,"coupler" => z_c)
        push!(config, "reconf" => z_l)
        push!(substation_configs, substation => config)
    end
    return solve_NCUC_L2!(grid,α_gt, β_gt,time_limit,T_=T_,reference_gen=reference_gen,warm_start=substation_configs)
end

function solve_NCUC_L2_SP(grid ::PowerGrid,α_gt,β_gt,substation,time_limit; T_ = nothing, reference_gen=nothing)
    ## NOT FINISHED YET
    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_

    Nodes_set = keys(grid.Buses)
    aux_bus_set = grid.Substations[substation].Aux_Buses_IDs

    my_coupler_id = first(grid.Substations[substation].Reconf_CouplerLines_IDs)
    my_reconf_ids = grid.Substations[substation].Reconf_AuxLines_IDs

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]

    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    other_couplers = collect(setdiff(Set(Coupler_set),Set(my_coupler_id)))
    other_reconfs = collect(setdiff(Set(Reconf_set),Set(my_reconf_ids)))

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    normally_opened_reconf_lines_per_coupler = Dict()

    for substation_id in keys(grid.Substations)
        substation_obj = grid.Substations[substation_id]
        coupler_id = substation_obj.Reconf_CouplerLines_IDs[1]
        reconf_ids = substation_obj.Reconf_AuxLines_IDs
        normally_opened_ids = []
        for reconf_id in reconf_ids
            if grid.Branches[reconf_id].GeneralSwitch.SwitchingStatus == 0
                push!(normally_opened_ids)
            end
        end
        push!(normally_opened_reconf_lines_per_coupler, coupler_id => normally_opened_ids)
    end

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", time_limit*0.5)
    # set_optimizer_attribute(m, "MIPGap", 0.01)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set, t in T] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] ≤ p[g in Gen_set,t in T] ≤ grid.Generators[g].Pg_max*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] )
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, t in T; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, z_c[coupler = Coupler_set, t in T], Bin)
    JuMP.@variable(m, z_l[reconf = Reconf_set, t in T], Bin)

    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    # 2. Constraints

    for t in T
        for coupler in my_coupler_id
            set_start_value(z_c[coupler, t], grid.Branches[coupler].GeneralSwitch.SwitchingStatus_t[t])
        end
        for reconf in my_reconf_ids
            set_start_value(z_l[reconf, t ], grid.Branches[reconf].GeneralSwitch.SwitchingStatus_t[t])
        end

        for coupler in other_couplers
            JuMP.fix(z_c[coupler,t],grid.Branches[coupler].GeneralSwitch.SwitchingStatus_t[t])
        end

        for reconf in other_reconfs
            JuMP.fix(z_l[reconf,t],grid.Branches[reconf].GeneralSwitch.SwitchingStatus_t[t])
        end
    end

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[t in T], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,t] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set, t in T],
        sum(pij[i,j,t] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        pij[i,j,t] == Sbase*(B[i,j])*(δ[i,t]-δ[j,t]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Branch_nodes], 
        pij[i,j,t] == -pij[j,i,t])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,t] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set, t in T],sum(z_l[l,t] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) ≤ 1)
    
    # 2.5.1 Normally switched of reconf lines if no splitting (lazy constraint)
    M = 1e3 # BigM number
    JuMP.@constraint(m,switched_off_reconf[coupler in Coupler_set, t in T], sum(z_l[l,t] for l in normally_opened_reconf_lines_per_coupler[coupler]) ≤ (1-z_c[coupler, t])*M)
    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[j,t]-(1-z_l[Reconf_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_l[Reconf_dict[ Set([i,j])],t])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[j,t] - (1-z_c[Coupler_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_c[Coupler_dict[ Set([i,j])],t]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        -z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        pij[i,j,t] ≤ z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        pij[i,j,t] ≤ z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))
    
    for c in my_coupler_id
        for t in T
            MOI.set(m, Gurobi.VariableAttribute("BranchPriority"), z_c[c,t], 1) 
        end
    end

    JuMP.optimize!(m)

    solution_status = raw_status(m)
    z_l_sol = Dict()
    for reconf in my_reconf_ids
        reconf_sol = Dict()
        for t in T
            push!(reconf_sol, t=>JuMP.value.(z_l[reconf,t]))
        end
        push!(z_l_sol,reconf => reconf_sol)
    end
    z_c_sol = Dict()
    for coupler in my_coupler_id
        coupler_sol = Dict()
        for t in T
            push!(coupler_sol, t => JuMP.value.(z_c[coupler,:]))
        end
        push!(z_c_sol,coupler => coupler_sol)
    end
    return z_l_sol, z_c_sol
end

function solve_NCUC_L1!(grid ::PowerGrid;T_=nothing,reference_gen=nothing)
    
    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]

    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Gurobi.Optimizer)
    # set_optimizer_attribute(m, "TimeLimit", 600)
    # set_optimizer_attribute(m, "MIPGap", 0.01)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set, t in T] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, p[g in Gen_set,t in T])
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, t in T; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, u_gt[g in Gen_set, t in T], Bin)
    JuMP.@variable(m, α_gt[g in Gen_set, t in T], Bin) # start-up
    JuMP.@variable(m, β_gt[g in Gen_set, t in T], Bin) # shut-down
    JuMP.@variable(m,z_c[coupler in Coupler_set,t in T], Bin)

    JuMP.@constraint(m, Pg_max[g in Gen_set, t in T], p[g,t] ≤ grid.Generators[g].Pg_max*u_gt[g,t])
    JuMP.@constraint(m, Pg_min[g in Gen_set, t in T], p[g,t] ≥ grid.Generators[g].Pg_min*u_gt[g,t])
    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    JuMP.@constraint(m, logic[g in Gen_set, t in T[2:T_]], u_gt[g,t]-u_gt[g,t-1] == α_gt[g,t] - β_gt[g,t])
    JuMP.@constraint(m, logic_init[g in Gen_set], u_gt[g,1] == α_gt[g,1] + β_gt[g,1])
    
    JuMP.@constraint(m, MUT[g in Gen_set, t in grid.Generators[g].min_up_time:T_],
        sum(α_gt[g,i] for i in t-grid.Generators[g].min_up_time+1:t) ≤ u_gt[g,t])

    JuMP.@constraint(m, MDT[g in Gen_set, t in grid.Generators[g].min_down_time:T_],
        sum(β_gt[g,i] for i in t-grid.Generators[g].min_down_time+1:t) ≤ 1-u_gt[g,t])
    
    # 2. Constraints

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[t in T], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,t] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set, t in T],
        sum(pij[i,j,t] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        pij[i,j,t] == Sbase*(B[i,j])*(δ[i,t]-δ[j,t]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Branch_nodes], 
        pij[i,j,t] == -pij[j,i,t])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,t] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)
    

    ######################

    JuMP.@constraint(m, s_phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[j,t]-(1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus_t[t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, s_phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[i,t] ≤ δ[j,t] + (1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus_t[t])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, s_phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[j,t] - (1-z_c[Coupler_dict[Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, s_phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_c[Coupler_dict[Set([i,j])],t]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,s_reconf_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        -grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus_t[t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,s_reconf_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        pij[i,j,t] ≤ grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus_t[t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,s_coupler_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,s_coupler_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        pij[i,j,t] ≤  z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*u_gt[g,t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))
    

    JuMP.optimize!(m)

    solution_status = raw_status(m)

    if has_values(m)
        grid.Operating_Cost = JuMP.objective_value(m)
        Pgt = JuMP.value.(p)
        Pijt = JuMP.value.(pij)
        δt = JuMP.value.(δ)
        u_gt = JuMP.value.(u_gt)

        for g in keys(grid.Generators)
            grid.Generators[g].GeneralSwitch.SwitchingStatus_t = [u_gt[g,t] for t in T]
            grid.Generators[g].Pg_t = [Pgt[g,t] for t in T]
            grid.Generators[g].Qg_t = zeros(1,T_)
        end

        for b in keys(grid.Branches)
            Fr_bus = grid.Branches[b].Fr_bus_ID
            To_bus = grid.Branches[b].To_bus_ID
            grid.Branches[b].PowerFlow_ij_t = [Pijt[Fr_bus,To_bus,t] for t in T]
            grid.Branches[b].PowerFlow_ji_t = [Pijt[To_bus,Fr_bus,t] for t in T]
            grid.Branches[b].losses_P_t = [Pijt[Fr_bus,To_bus,t]+Pijt[To_bus,Fr_bus,t] for t in T]
            
            grid.Branches[b].ReactFlow_ij_t = zeros(1,T_)
            grid.Branches[b].ReactFlow_ji_t = zeros(1,T_)
            grid.Branches[b].losses_Q_t = zeros(1,T_)

            branch_type = grid.Branches[b].BranchType
            if branch_type == 2 # coupler line
                grid.Branches[b].GeneralSwitch.SwitchingStatus_t = [JuMP.value.(z_c[b,t]) for t in T]
            end

        end

        for bus in keys(grid.Buses)
            grid.Buses[bus].δ_t = [δt[bus,t] for t in T]
            grid.Buses[bus].V_magnitude_t = ones(1,T_)
        end

        
        return JuMP.objective_value(m) , JuMP.value.(α_gt), JuMP.value.(β_gt)
    else
        return solution_status
    end
end

function solve_NCUC_L2!(grid ::PowerGrid,α_gt,β_gt, time_limit;T_=nothing,reference_gen=nothing,warm_start=nothing)
    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]

    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    normally_opened_reconf_lines_per_coupler = Dict()

    for substation_id in keys(grid.Substations)
        substation_obj = grid.Substations[substation_id]
        coupler_id = substation_obj.Reconf_CouplerLines_IDs[1]
        reconf_ids = substation_obj.Reconf_AuxLines_IDs
        normally_opened_ids = []
        for reconf_id in reconf_ids
            if grid.Branches[reconf_id].GeneralSwitch.SwitchingStatus == 0
                push!(normally_opened_ids)
            end
        end
        push!(normally_opened_reconf_lines_per_coupler, coupler_id => normally_opened_ids)
    end

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", time_limit)
    # set_optimizer_attribute(m, "MIPGap", 0.01)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set, t in T] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] ≤ p[g in Gen_set,t in T] ≤ grid.Generators[g].Pg_max*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] )
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, t in T; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, z_c[coupler = Coupler_set, t in T], Bin)
    JuMP.@variable(m, z_l[reconf = Reconf_set, t in T], Bin)

    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    # 2. Constraints

    if warm_start === nothing
        for t in T
            for coupler in Coupler_set
                set_start_value(z_c[coupler, t], grid.Branches[coupler].GeneralSwitch.SwitchingStatus_t[t])
            end
            for reconf in Reconf_set
                set_start_value(z_l[reconf, t ], grid.Branches[reconf].GeneralSwitch.SwitchingStatus_t[t])
            end
        end
    else
        for substation_id in keys(warm_start)
            substation = warm_start[substation_id]
            my_coupler_id = first(grid.Substations[substation_id].Reconf_CouplerLines_IDs)
            my_reconf_ids = grid.Substations[substation_id].Reconf_AuxLines_IDs
            for t in T
                for coupler in my_coupler_id
                    set_start_value(z_c[coupler, t], first(substation["coupler"][coupler][t]))
                end
                for reconf in my_reconf_ids
                    set_start_value(z_l[reconf, t], first(substation["reconf"][reconf][t]))
                end
            end
        end
    end

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[t in T], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,t] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set, t in T],
        sum(pij[i,j,t] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        pij[i,j,t] == Sbase*(B[i,j])*(δ[i,t]-δ[j,t]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Branch_nodes], 
        pij[i,j,t] == -pij[j,i,t])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,t] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set, t in T], sum(z_l[l,t] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
    
    # 2.5.1 Normally switched off reconf lines if no splitting (lazy constraint)
    M = 1e3 # BigM number
    JuMP.@constraint(m,switched_off_reconf[coupler in Coupler_set, t in T], sum(z_l[l,t] for l in normally_opened_reconf_lines_per_coupler[coupler]) ≤ (1-z_c[coupler, t])*M)
    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[j,t]-(1-z_l[Reconf_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_l[Reconf_dict[ Set([i,j])],t])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[j,t] - (1-z_c[Coupler_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes],
        δ[i,t] ≤ δ[j,t] + (1-z_c[Coupler_dict[ Set([i,j])],t]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        -z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Reconf_nodes], 
        pij[i,j,t] ≤ z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set, t in T; Set([i,j]) in Coupler_nodes], 
        pij[i,j,t] ≤ z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))
    
    for c in Coupler_set
        for t in T
            MOI.set(m, Gurobi.VariableAttribute("BranchPriority"), z_c[c,t], 1) 
        end
    end

    JuMP.optimize!(m)

    solution_status = raw_status(m)

    if has_values(m)
        grid.Operating_Cost = JuMP.objective_value(m)
        Pgt = JuMP.value.(p)
        Pijt = JuMP.value.(pij)
        δt = JuMP.value.(δ)
        Z_c = JuMP.value.(z_c)
        Z_l = JuMP.value.(z_l)

        for g in keys(grid.Generators)
            grid.Generators[g].Pg_t = [Pgt[g,t] for t in T]
            grid.Generators[g].Qg_t = zeros(1,T_)
        end

        for b in keys(grid.Branches)
            Fr_bus = grid.Branches[b].Fr_bus_ID
            To_bus = grid.Branches[b].To_bus_ID
            grid.Branches[b].PowerFlow_ij_t = [Pijt[Fr_bus,To_bus,t] for t in T]
            grid.Branches[b].PowerFlow_ji_t = [Pijt[To_bus,Fr_bus,t] for t in T]
            grid.Branches[b].losses_P_t = [Pijt[Fr_bus,To_bus,t]+Pijt[To_bus,Fr_bus,t] for t in T]
            
            grid.Branches[b].ReactFlow_ij_t = zeros(1,T_)
            grid.Branches[b].ReactFlow_ji_t = zeros(1,T_)
            grid.Branches[b].losses_Q_t = zeros(1,T_)

            branch_type = grid.Branches[b].BranchType
            if branch_type == 1 # reconf line
                grid.Branches[b].GeneralSwitch.SwitchingStatus_t = [Z_l[b,t] for t in T]
            elseif branch_type == 2 # coupler line
                grid.Branches[b].GeneralSwitch.SwitchingStatus_t = [Z_c[b,t] for t in T]
            end
        end

        for bus in keys(grid.Buses)
            grid.Buses[bus].δ_t = [δt[bus,t] for t in T]
            grid.Buses[bus].V_magnitude_t = ones(1,T_)
        end

        return JuMP.objective_value(m)
    else
        return solution_status
    end
end

function solve_UC_MP!(grid ::PowerGrid,t_reconf_set;T_=nothing,reference_gen=nothing)
    
    if T_ === nothing
        T_ = grid.N_time_steps
    else
        T_ = minimum([T_,grid.N_time_steps])
    end

    T = 1:T_

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]

    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", 300)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set, t in T] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, p[g in Gen_set,t in T])
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, t in T; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, u_gt[g in Gen_set, t in T], Bin)
    JuMP.@variable(m, α_gt[g in Gen_set, t in T], Bin) # start-up
    JuMP.@variable(m, β_gt[g in Gen_set, t in T], Bin) # shut-down
    if ! isempty(collect(t_reconf_set))
        JuMP.@variable(m, z_c[coupler = Coupler_set, t in collect(t_reconf_set)], Bin)
        JuMP.@variable(m, z_l[reconf = Reconf_set, t in collect(t_reconf_set)], Bin)
        for t in collect(t_reconf_set)
            for coupler in Coupler_set
                set_start_value(z_c[coupler, t], grid.Branches[coupler].GeneralSwitch.SwitchingStatus_t[t])
            end
            for reconf in Reconf_set
                set_start_value(z_l[reconf, t ], grid.Branches[reconf].GeneralSwitch.SwitchingStatus_t[t])
            end
        end
        
        
    end

    JuMP.@constraint(m, Pg_max[g in Gen_set, t in T], p[g,t] ≤ grid.Generators[g].Pg_max*u_gt[g,t])
    JuMP.@constraint(m, Pg_min[g in Gen_set, t in T], p[g,t] ≥ grid.Generators[g].Pg_min*u_gt[g,t])
    JuMP.@constraint(m, Ramp_Up[g in Gen_set, t in T[2:T_]], p[g,t]-p[g,t-1] ≤ grid.Generators[g].Δ_up)
    JuMP.@constraint(m, Ramp_Down[g in Gen_set, t in T[2:T_]], p[g,t-1]-p[g,t] ≤ grid.Generators[g].Δ_down)

    JuMP.@constraint(m, logic[g in Gen_set, t in T[2:T_]], u_gt[g,t]-u_gt[g,t-1] == α_gt[g,t] - β_gt[g,t])
    JuMP.@constraint(m, logic_init[g in Gen_set], u_gt[g,1] == α_gt[g,1] + β_gt[g,1])
    
    JuMP.@constraint(m, MUT[g in Gen_set, t in grid.Generators[g].min_up_time:T_],
        sum(α_gt[g,i] for i in t-grid.Generators[g].min_up_time+1:t) ≤ u_gt[g,t])

    JuMP.@constraint(m, MDT[g in Gen_set, t in grid.Generators[g].min_down_time:T_],
        sum(β_gt[g,i] for i in t-grid.Generators[g].min_down_time+1:t) ≤ 1-u_gt[g,t])
    
    # 2. Constraints

    # JuMP.@constraint(m, initial_reconf_state[l in Reconf_set,t in collect(setdiff(Set(T),t_reconf_set))],
    #     z_l[l,t] == grid.Branches[l].GeneralSwitch.SwitchingStatus)
    
    # JuMP.@constraint(m, initial_split_state[c in Coupler_set,t in collect(setdiff(Set(T),t_reconf_set))],
    #     z_c[c,t] == grid.Branches[c].GeneralSwitch.SwitchingStatus)

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[t in T], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,t] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set, t in T],
        sum(pij[i,j,t] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        pij[i,j,t] == Sbase*(B[i,j])*(δ[i,t]-δ[j,t]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Branch_nodes], 
        pij[i,j,t] == -pij[j,i,t])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set, t in T; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,t] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)
    
    if ! isempty(collect(t_reconf_set))
        # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
        JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set, t in collect(t_reconf_set)],sum(z_l[l,t] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)

        # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Reconf_nodes],
            δ[j,t]-(1-z_l[Reconf_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

        JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Reconf_nodes],
            δ[i,t] ≤ δ[j,t] + (1-z_l[Reconf_dict[ Set([i,j])],t])) 

        # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Coupler_nodes],
            δ[j,t] - (1-z_c[Coupler_dict[ Set([i,j])],t]) ≤ δ[i,t] ) 

        JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Coupler_nodes],
            δ[i,t] ≤ δ[j,t] + (1-z_c[Coupler_dict[ Set([i,j])],t]) ) 

        # 2.7.1 Reconfiguration line capacity
        JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Reconf_nodes], 
            -z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

        JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Reconf_nodes], 
            pij[i,j,t] ≤ z_l[Reconf_dict[Set([i,j])],t] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

        # 2.7.2 Reconfiguration line capacity
        JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Coupler_nodes], 
            -z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

        JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set, t in collect(t_reconf_set); Set([i,j]) in Coupler_nodes], 
            pij[i,j,t] ≤ z_c[Coupler_dict[Set([i,j])],t] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)
    end

    ######################

    JuMP.@constraint(m, s_phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Reconf_nodes],
        δ[j,t]-(1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ≤ δ[i,t] ) 

    JuMP.@constraint(m, s_phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Reconf_nodes],
        δ[i,t] ≤ δ[j,t] + (1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus)) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, s_phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Coupler_nodes],
        δ[j,t] - (1-grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ≤ δ[i,t] ) 

    JuMP.@constraint(m, s_phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Coupler_nodes],
        δ[i,t] ≤ δ[j,t] + (1-grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,s_reconf_cap_1[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Reconf_nodes], 
        -grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,s_reconf_cap_2[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Reconf_nodes], 
        pij[i,j,t] ≤ grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,s_coupler_cap_1[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Coupler_nodes], 
        -grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,t])

    JuMP.@constraint(m,s_coupler_cap_2[i in Nodes_set, j in Nodes_set, t in collect(setdiff(Set(T),t_reconf_set)); Set([i,j]) in Coupler_nodes], 
        pij[i,j,t] ≤  grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,t]+grid.Generators[g].C0*u_gt[g,t]+α_gt[g,t]*grid.Generators[g].start_up_cost+β_gt[g,t]*grid.Generators[g].shut_down_cost for g in Gen_set, t in T))

    JuMP.optimize!(m)

    solution_status = raw_status(m)

    if has_values(m)
        grid.Operating_Cost = JuMP.objective_value(m)
        Pgt = JuMP.value.(p)
        Pijt = JuMP.value.(pij)
        δt = JuMP.value.(δ)
        u_gt = JuMP.value.(u_gt)

        for g in keys(grid.Generators)
            grid.Generators[g].GeneralSwitch.SwitchingStatus_t = [u_gt[g,t] for t in T]
            grid.Generators[g].Pg_t = [Pgt[g,t] for t in T]
            grid.Generators[g].Qg_t = zeros(1,T_)
        end

        for b in keys(grid.Branches)
            Fr_bus = grid.Branches[b].Fr_bus_ID
            To_bus = grid.Branches[b].To_bus_ID
            grid.Branches[b].PowerFlow_ij_t = [Pijt[Fr_bus,To_bus,t] for t in T]
            grid.Branches[b].PowerFlow_ji_t = [Pijt[To_bus,Fr_bus,t] for t in T]
            grid.Branches[b].losses_P_t = [Pijt[Fr_bus,To_bus,t]+Pijt[To_bus,Fr_bus,t] for t in T]
            
            grid.Branches[b].ReactFlow_ij_t = zeros(1,T_)
            grid.Branches[b].ReactFlow_ji_t = zeros(1,T_)
            grid.Branches[b].losses_Q_t = zeros(1,T_)

            branch_type = grid.Branches[b].BranchType
            if branch_type == 1 # reconf line
                for t in collect(t_reconf_set)
                    grid.Branches[b].GeneralSwitch.SwitchingStatus_t[t] = JuMP.value.(z_l[b,t])
                end
            elseif branch_type == 2 # coupler line
                for t in collect(t_reconf_set)
                    grid.Branches[b].GeneralSwitch.SwitchingStatus_t[t] = JuMP.value.(z_c[b,t])
                end
            end
        end

        for bus in keys(grid.Buses)
            grid.Buses[bus].δ_t = [δt[bus,t] for t in T]
            grid.Buses[bus].V_magnitude_t = ones(1,T_)
        end

        t_l_congested_dict = Dict()
        for t in T
            for l in Transmission_set
                Pijt = abs(grid.Branches[l].PowerFlow_ij_t[t])
                rate_l = grid.Branches[l].rating*grid.S_base
                if Pijt ≥ rate_l
                    if t ∉ keys(t_l_congested_dict)
                        push!(t_l_congested_dict, t => [l])
                    else
                        push!(t_l_congested_dict[t],l)
                    end
                end
            end
        end

        T_congested = Set([t for t in T if any(abs(grid.Branches[l].PowerFlow_ij_t[t])==grid.Branches[l].rating*grid.S_base for l in Transmission_set)])
        @assert T_congested == Set(collect(keys(t_l_congested_dict)))
        return JuMP.objective_value(m) , T_congested, t_l_congested_dict
    end
end

function solve_UC_SP_OBS!(grid ::PowerGrid, T_congested_new, UB; reference_gen=nothing)

    reduced_cost = Dict()
    for t in collect(T_congested_new)
        UB_t = sum(grid.Generators[g].Pg_t[t]*grid.Generators[g].C1 + grid.Generators[g].C0 for g in keys(grid.Generators) if grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] == 1)
        cost = solve_OBS_t!(grid, t, reference_gen=reference_gen)
        Δ = UB_t - cost
        push!(reduced_cost,Δ => t)
    end
    max_reduction = maximum(collect(keys(reduced_cost)))
    t_reconf = reduced_cost[max_reduction]
    LB = UB - max_reduction
    return LB, t_reconf
end

function solve_UC_SP_MC!(grid ::PowerGrid, T_congested_new, t_l_congested_dict, UB; reference_gen=nothing)

    total_load_shift = Dict()
    for t in collect(T_congested_new)
        load_shift = solve_min_congestion_t!(grid, t, t_l_congested_dict[t], reference_gen=reference_gen)
        push!(total_load_shift,load_shift => t)
    end
    max_load_shift = maximum(collect(keys(total_load_shift)))
    t_reconf = total_load_shift[max_load_shift]
    UB_t = sum(grid.Generators[g].Pg_t[t_reconf]*grid.Generators[g].C1 + grid.Generators[g].C0 for g in keys(grid.Generators) if grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t_reconf] == 1)
    max_reduction = UB_t - solve_ED_t!(grid,t_reconf)
    LB = UB - max_reduction
    return LB, t_reconf
end

function solve_UC_SP_ED!(grid ::PowerGrid, T_congested_new, UB)

    reduced_cost = Dict()
    for t in collect(T_congested_new)
        UB_t = sum(grid.Generators[g].Pg_t[t]*grid.Generators[g].C1 + grid.Generators[g].C0 for g in keys(grid.Generators) if grid.Generators[g].GeneralSwitch.SwitchingStatus_t[t] == 1)
        cost = solve_ED_t!(grid, t)
        Δ = UB_t - cost
        push!(reduced_cost,Δ => t)
    end
    max_reduction = maximum(collect(keys(reduced_cost)))
    t_reconf = reduced_cost[max_reduction]
    LB = UB - max_reduction
    return LB, t_reconf
end

function solve_OBS_t!(grid ::PowerGrid, t; reference_gen=nothing)
    
    m = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, "TimeLimit", 60)

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = [k for k in keys(grid.Generators) if grid.Generators[k].GeneralSwitch.SwitchingStatus_t[t] == 1]
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, z_c[coupler = Coupler_set], Bin)
    JuMP.@variable(m, z_l[reconf = Reconf_set], Bin)

    # 2. Constraints

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in Transmission_nodes],
        pij[i,j] == Sbase*(B[i,j])*(δ[i]-δ[j]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes], 
        pij[i,j] == -pij[j,i])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set],sum(z_l[l] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) ≤ 1)

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[j]-(1-z_l[Reconf_dict[ Set([i,j])]]) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[i] ≤ δ[j] + (1-z_l[Reconf_dict[ Set([i,j])]])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[j] - (1-z_c[Coupler_dict[ Set([i,j])]]) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[i] ≤ δ[j] + (1-z_c[Coupler_dict[ Set([i,j])]]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        -z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        pij[i,j] ≤ z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        pij[i,j] ≤ z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))

    JuMP.optimize!(m)

    # for b in keys(grid.Branches)

    #     branch_type = grid.Branches[b].BranchType
    #     if branch_type == 1 # reconf line
    #         grid.Branches[b].GeneralSwitch.SwitchingStatus_t[t] = JuMP.value.(z_l[b])
    #     elseif branch_type == 2 # coupler line
    #         grid.Branches[b].GeneralSwitch.SwitchingStatus_t[t] = JuMP.value.(z_c[b])
    #     end
    # end

    return JuMP.objective_value(m)
end

function solve_min_congestion_t!(grid ::PowerGrid,t,l_congested_arr; reference_gen=nothing)

    Nodes_set = keys(grid.Buses)
    aux_bus_set = [key for key in Nodes_set if grid.Buses[key].BusType == 1]

    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_set = [key for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_set = [key for key in Branch_set if grid.Branches[key].BranchType == 2]

    Reconf_dict = Dict()
    for key in Reconf_set
        push!(Reconf_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Coupler_dict = Dict()
    for key in Coupler_set
        push!(Coupler_dict, Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) => key)
    end

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 0]
    Reconf_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 1]
    Coupler_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Branch_set if grid.Branches[key].BranchType == 2]

    Gen_set = [k for k in keys(grid.Generators) if grid.Generators[k].GeneralSwitch.SwitchingStatus_t[t] == 1]
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    m = JuMP.Model(Ipopt.Optimizer)
    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, LS[d in keys(grid.Loads)])

    # 2. Constraints

    # JuMP.@constraint(m,load_shedding_limits[d in keys(grid.Loads)], -grid.Loads[d].Pd_t[t] ≤ LS[d] ≤ grid.Loads[d].Pd_t[t])
    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    JuMP.@constraint(m, Nodal_balance[i in Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(grid.Generators[g].Pg_t[t] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i) + sum(LS[d] for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

    JuMP.@constraint(m, net_zero_shedding, sum(LS[d] for d in keys(grid.Loads)) == 0)
    # 2.3.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in Transmission_nodes],
        pij[i,j] == Sbase*(B[i,j])*(δ[i]-δ[j]))

    # 2.3.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes], 
        pij[i,j] == -pij[j,i])

    # 2.4 Thermal limits of transmission lines
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[j]-(1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[i] ≤ δ[j] + (1-grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus)) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[j] - (1-grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[i] ≤ δ[j] + (1-grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        -grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        pij[i,j] ≤ grid.Branches[Reconf_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        -grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        pij[i,j] ≤  grid.Branches[Coupler_dict[Set([i,j])]].GeneralSwitch.SwitchingStatus * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(pij[grid.Branches[l].Fr_bus_ID,grid.Branches[l].To_bus_ID]^2 for l in l_congested_arr))
    # JuMP.@objective(m,Min,sum(pij[collect(l)]^2-Sbase*grid.Branches[l].rating^2 for l in Transmission_nodes))

    JuMP.optimize!(m)

    total_load_shifting = sum(JuMP.value.(LS[d]) for d in keys(grid.Loads) if JuMP.value.(LS[d]) ≥ 0)
    return total_load_shifting
end

function solve_ED_t!(grid ::PowerGrid,t)
    m = JuMP.Model(Gurobi.Optimizer)
    Gen_set = [k for k in keys(grid.Generators) if grid.Generators[k].GeneralSwitch.SwitchingStatus_t[t] == 1]
    # 1. Variables
    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    # 2. Constraints
    JuMP.@constraint(m,Nodal_balance,
        sum(p[g] for g in Gen_set) == sum(grid.Loads[d].Pd_t[t] for d in keys(grid.Loads)))

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))

    JuMP.optimize!(m)

    return JuMP.objective_value(m)
end