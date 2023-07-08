function solve_OPF!(grid ::PowerGrid,method, reference_gen=nothing)

    if method == :ACOPF
        return _solve_ACOPF!(grid,reference_gen)
    elseif method == :DCOPF
        return _solve_DCOPF!(grid,reference_gen)
    end
    
end

function _solve_DCOPF!(grid ::PowerGrid,reference_gen=nothing)
    m = JuMP.Model(Gurobi.Optimizer)

    Nodes_set = keys(grid.Buses)
    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end
    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])

    # 2. Constraints

    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Pd)

    JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
        pij[i,j] == Sbase*(B[i,j])*(δ[i]-δ[j]))

    JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))

    JuMP.optimize!(m)
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."
        Line_Duals = Dict()
        for branch in Branch_set
            current_dual1 = JuMP.dual(pl_rate[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID])
            current_dual2 = JuMP.dual(pl_rate[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID])
            Line_Duals = push!(Line_Duals,branch => (current_dual1,current_dual2))
        end

        Bus_Duals = Dict()
        for bus in Nodes_set
            current_dual_b = JuMP.dual(Nodal_balance[bus])
            Bus_Duals = push!(Bus_Duals,bus => current_dual_b)
        end

        grid.Line_Duals = Line_Duals
        grid.Bus_Duals = Bus_Duals

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)

        [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
        [grid.Generators[g].Qg = 0 for g in Gen_set]

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        update_grid_tables!(grid)
    end

    return solution_status
    
end

function _solve_ACOPF!(grid ::PowerGrid,reference_gen=nothing)
    m = JuMP.Model(Ipopt.Optimizer)

    Nodes_set = keys(grid.Buses)
    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end
    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    B = imag(grid.Y_bus)
    G = real(grid.Y_bus)


    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].V_min ≤ v[i in Nodes_set] ≤ grid.Buses[i].V_max)
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)

    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    JuMP.@variable(m, grid.Generators[g].Qg_min ≤ q[g in Gen_set] ≤ grid.Generators[g].Qg_max)

    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, qij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])

    # 2. Constraints
    
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    # ACTIVE POWER THROUGH LINE i-j
    JuMP.@NLconstraint(m, p_line[i in Nodes_set, j in Nodes_set; Set([i,j]) in Branch_nodes],
    pij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*cos(δ[i]-δ[j])+B[i,j]*sin(δ[i]-δ[j])) -(v[i]^2)*G[i,j]))

    # REACTIVE POWER THROUGH LINE i-j
    JuMP.@NLconstraint(m, q_line[i in Nodes_set, j in Nodes_set; Set([i,j]) in Branch_nodes],
    qij[i,j] ==  Sbase*(v[i]*v[j]*(G[i,j]*sin(δ[i]-δ[j]) - B[i,j]*cos(δ[i]-δ[j])) +(v[i]^2)*(B[i,j] - grid.Branches[branch_dictionary[Set([i,j])]].b/2)))

    # ACTIVE NODAL BALANCE
    JuMP.@constraint(m, Pnodal[i in Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Pd)

    # REACTIVE NODAL BALANCE
    JuMP.@constraint(m, Qnodal[i in Nodes_set],
        sum(qij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(q[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Qd)

    # LINE CAPACITY
    JuMP.@NLconstraint(m, Smax[i in Nodes_set, j in Nodes_set; Set([i,j]) in Branch_nodes],
        pij[i,j]^2 + qij[i,j]^2 ≤ ((Sbase)^2)*(grid.Branches[branch_dictionary[Set([i,j])]].rating)^2)

    # 3. OBJECTIVE
    JuMP.@NLobjective(m,Min,sum(grid.Generators[g].C2*p[g]^2+grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))

    JuMP.optimize!(m)
    solution_status = raw_status(m)

    if solution_status != "Infeasible_Problem_Detected"

        Line_Duals = Dict()
        for branch in Branch_set
            current_dual1 = JuMP.dual(Smax[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID])
            current_dual2 = JuMP.dual(Smax[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID])
            Line_Duals = push!(Line_Duals,branch => (current_dual1,current_dual2))
        end

        Bus_Duals = Dict()
        for bus in Nodes_set
            current_dual_b = JuMP.dual(Pnodal[bus])
            Bus_Duals = push!(Bus_Duals,bus => current_dual_b)
        end

        grid.Line_Duals = Line_Duals
        grid.Bus_Duals = Bus_Duals

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)

        [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
        [grid.Generators[g].Qg = JuMP.value.(q[g]) for g in Gen_set]

        [grid.Buses[bus].V_magnitude = JuMP.value.(v[bus]) for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = JuMP.value.(qij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = JuMP.value.(qij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_Q = grid.Branches[branch].ReactFlow_ij + grid.Branches[branch].ReactFlow_ji for branch in Branch_set]

        update_grid_tables!(grid)
    end

    return solution_status

end

function solve_OTS!(grid ::PowerGrid,reference_gen=nothing, max_op = Inf)
    m = JuMP.Model(Gurobi.Optimizer)

    Nodes_set = keys(grid.Buses)
    Branch_set = keys(grid.Branches)
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end
    Gen_set = keys(grid.Generators)
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m,z[branch_id in keys(grid.Branches)], Bin)

    # 2. Constraints

    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    JuMP.@constraint(m,Nodal_balance[i in Nodes_set],
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Pd)
    
    M = Sbase*100
    JuMP.@constraint(m,pl_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
        pij[i,j] - Sbase*(B[i,j])*(δ[i]-δ[j]) ≤  (1-z[branch_dictionary[Set([i,j])]])*M)
    
    JuMP.@constraint(m,pl_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
        pij[i,j] - Sbase*(B[i,j])*(δ[i]-δ[j]) ≥  -(1-z[branch_dictionary[Set([i,j])]])*M)

    JuMP.@constraint(m,pl_rate_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])]] ≤ pij[i,j] )

    JuMP.@constraint(m,pl_rate_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes],
         pij[i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])]] )


    if max_op !== Inf
        JuMP.@constraint(m,maximum_allowed_switching_actions, sum(1-z[l] for l in Branch_set) ≤ max_op)
    end

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))

    JuMP.optimize!(m)
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)

        [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
        [grid.Generators[g].Qg = 0 for g in Gen_set]

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        [grid.Branches[branch].GeneralSwitch.SwitchingStatus = JuMP.value.(z[branch]) for branch in Branch_set]

        grid.Z_lines = JuMP.value.(z)
        reduce_grid!(grid)
        update_grid_tables!(grid)
    end

    return solution_status
end

function solve_OBS!(grid ::PowerGrid,reference_gen=nothing)

    m = JuMP.Model(Gurobi.Optimizer)

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
        sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

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

    M_δ = 2*π
    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[j]-M_δ*(1-z_l[Reconf_dict[ Set([i,j])]]) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[i] ≤ δ[j] + M_δ*(1-z_l[Reconf_dict[ Set([i,j])]])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[j] - M_δ*(1-z_c[Coupler_dict[ Set([i,j])]]) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[i] ≤ δ[j] + M_δ*(1-z_c[Coupler_dict[ Set([i,j])]]) ) 

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
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)

        [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
        [grid.Generators[g].Qg = 0 for g in Gen_set]

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        [grid.Branches[Reconf_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_l[Reconf_id]) for Reconf_id in Reconf_set]
        [grid.Branches[coupler_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_c[coupler_id]) for coupler_id in Coupler_set]
        
        grid.Z_lines = JuMP.value.(z_l)
        grid.Z_coupler = JuMP.value.(z_c)

        for bus in aux_bus_set
            if sum(JuMP.value.(z_l[l]) for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) < 1
                grid.Buses[bus].GeneralSwitch.SwitchingStatus = 0
            end
        end

        reduce_grid!(grid)
        update_grid_tables!(grid)
    end

    return solution_status
    
end

function solve_HOBS!(grid ::PowerGrid,reference_gen=nothing)
    m = JuMP.Model(Gurobi.Optimizer)

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

    virtual_generator_couplers = Dict()

    for key in keys(Coupler_dict)
        coupler_buses = collect(key)
        for bus in coupler_buses
            if grid.Buses[bus].ConnectedGensIDs != []
                push!(virtual_generator_couplers, grid.Buses[bus].ConnectedGensIDs[1] => Coupler_dict[key]) 
            end
        end
    end
  
    Gen_set = [key for key in keys(grid.Generators) if grid.Generators[key].GenType != :virtual]
    Virtual_Gen_set = [key for key in keys(grid.Generators) if grid.Generators[key].GenType == :virtual]
    Sbase = grid.S_base

    B = -1*imag(grid.Y_bus)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    if Virtual_Gen_set != []
        JuMP.@variable(m, p_vsc[g in Virtual_Gen_set])
    end
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set; Set([i,j]) in Branch_nodes])
    JuMP.@variable(m, z_c[coupler = Coupler_set], Bin)
    JuMP.@variable(m, z_l[reconf = Reconf_set], Bin)

    # 2. Constraints

    # 2.0 virtual generator settings
    if Virtual_Gen_set !=[]
        JuMP.@constraint(m,VSC_cap1[g in Virtual_Gen_set],(1-z_c[virtual_generator_couplers[g]])*grid.Generators[g].Pg_min ≤ p_vsc[g])
        JuMP.@constraint(m,VSC_cap2[g in Virtual_Gen_set], p_vsc[g] ≤ (1-z_c[virtual_generator_couplers[g]])*grid.Generators[g].Pg_max)
        JuMP.@constraint(m,VSC[l in Coupler_set],p_vsc[grid.Buses[grid.Branches[l].Fr_bus_ID].ConnectedGensIDs[1]] == -p_vsc[grid.Buses[grid.Branches[l].To_bus_ID].ConnectedGensIDs[1]])
    end

    # 2.1 Refernece angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle, δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID] ==  0.0)
    end

    # 2.2 Nodal balance for all nodes
    if Virtual_Gen_set != []
        JuMP.@constraint(m, Nodal_balance[i in Nodes_set],
            sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p_vsc[g] for g in Virtual_Gen_set if grid.Generators[g].GenBus_ID == i) + sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
    else
        JuMP.@constraint(m, Nodal_balance[i in Nodes_set],
            sum(pij[i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
    end
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
    JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set],sum(z_l[l] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)

    M_δ = 2*π
    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[j]-M_δ*(1-z_l[Reconf_dict[ Set([i,j])]]) ≤ δ[i]) 

    JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        δ[i] ≤ δ[j] + M_δ*(1-z_l[Reconf_dict[ Set([i,j])]])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[j] - M_δ*(1-z_c[Coupler_dict[ Set([i,j])]]) ≤ δ[i] ) 

    JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        δ[i] ≤ δ[j] + M_δ*(1-z_c[Coupler_dict[ Set([i,j])]]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        -z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        pij[i,j] ≤ z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Coupler line capacity
    JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        -z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j])

    JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        pij[i,j] ≤ z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    # 3. Objective
    JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set)-sum(z_c[c] for c in Coupler_set)*1e-3)

    JuMP.optimize!(m)
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)

        [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
        [grid.Generators[g].Qg = 0 for g in Gen_set]

        if Virtual_Gen_set != []
            [grid.Generators[g].Pg = JuMP.value.(p_vsc[g]) for g in Virtual_Gen_set]
            [grid.Generators[g].Qg = 0 for g in Virtual_Gen_set]
        end

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        [grid.Branches[Reconf_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_l[Reconf_id]) for Reconf_id in Reconf_set]
        [grid.Branches[coupler_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_c[coupler_id]) for coupler_id in Coupler_set]
        
        grid.z_reconf = JuMP.value.(z_l)
        grid.z_coupler = JuMP.value.(z_c)

        for bus in aux_bus_set
            if sum(JuMP.value.(z_l[l]) for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) < 1
                grid.Buses[bus].GeneralSwitch.SwitchingStatus = 0
            end
        end

        reduce_grid!(grid)
        update_grid_tables!(grid)
    end

    return solution_status
end
