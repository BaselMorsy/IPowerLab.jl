function solve_SCOPF_TC!(grid ::PowerGrid,Redispatch=false,method=:DCOPF,reference_gen=nothing)

    _,a,k = Generate_Contingencies!(grid,:TC)
    
    if method == :DCOPF
        return _solve_DC_SCOPF_TC!(grid,Redispatch,reference_gen,a,k)
    elseif method == :ACOPF
        return _solve_AC_SCOPF_TC!(grid,Redispatch,reference_gen,a,k)
    else
        println("This method has not been implemented yet!")
    end
end

function solve_SCOPF_TGC!(grid ::PowerGrid,method=:DCOPF)
    if method == :DCOPF
        return _solve_DC_SCOPF_TGC!(grid)
    elseif method == :ACOPF
        return _solve_AC_SCOPF_TGC!(grid)
    else
        println("This method has not been implemented yet!")
    end
end

function _solve_DC_SCOPF_TC!(grid ::PowerGrid,Redispatch=false,reference_gen=nothing,a=[],c=[])
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
    k_set = collect(1:c)

    B = -1*imag(grid.Y_bus)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set,k in k_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, k in k_set; Set([i,j]) in Branch_nodes])
    if Redispatch
        JuMP.@variable(m, grid.Generators[g,k].Pg_min ≤ p[g in Gen_set, k in k_set] ≤ grid.Generators[g,k].Pg_max)
        JuMP.@constraint(m,ramp_up[g in Gen_set,k in k_set[2:c]], p[g,k]-p[g,0] ≤ grid.Generators[g].Δ_up)
        JuMP.@constraint(m,ramp_up[g in Gen_set,k in k_set[2:c]], p[g,0]-p[g,k] ≤ grid.Generators[g].Δ_down)
    else
        JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    end
    

    # 2. Constraints

    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[k in k_set], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,k] ==  0.0)
    end

    if Redispatch
        JuMP.@constraint(m,Nodal_balance[i in Nodes_set,k in k_set],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Pd)
    else
        JuMP.@constraint(m,Nodal_balance[i in Nodes_set,k in k_set],
            sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - grid.Loads[i].Pd)
    end

    JuMP.@constraint(m,pl[i in Nodes_set,j in Nodes_set,k in k_set; Set([i,j]) in Branch_nodes],
        pij[i,j,k] == Sbase*(B[i,j])*(δ[i,k]-δ[j,k]))

    JuMP.@constraint(m,pl_rate[i in Nodes_set,j in Nodes_set,k in k_set; Set([i,j]) in Branch_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*a[i,j,k] ≤ pij[i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*a[i,j,k])

    # 3. Objective
    if Redispatch
        JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,1]+grid.Generators[g].C0 for g in Gen_set))
    else
        JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))
    end

    JuMP.optimize!(m)
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."
        Line_Duals = Dict()
        for branch in Branch_set
            current_dual1 = JuMP.dual(pl_rate[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,!])
            current_dual2 = JuMP.dual(pl_rate[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,!])
            Line_Duals = push!(Line_Duals,branch => (current_dual1,current_dual2))
        end

        Bus_Duals = Dict()
        for bus in Nodes_set
            current_dual_b = JuMP.dual(Nodal_balance[bus,!])
            Bus_Duals = push!(Bus_Duals,bus => current_dual_b)
        end

        grid.Line_Duals = Line_Duals
        grid.Bus_Duals = Bus_Duals

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)
        if Redispatch
            [grid.Generators[g].Pg = JuMP.value.(p[g,0]) for g in Gen_set]
            [grid.Generators[g].Qg = 0 for g in Gen_set]
            [grid.Generators[g].Pg_k = JuMP.value.(p[g,!]) for g in Gen_set]
            [grid.Generators[g].Qg_k = 0 for g in Gen_set]
        else
            [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
            [grid.Generators[g].Qg = 0 for g in Gen_set]
            [grid.Generators[g].Pg_k = JuMP.value.(p[g]) for g in Gen_set]
            [grid.Generators[g].Qg_k = 0 for g in Gen_set]
        end

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus,0]) for bus in Nodes_set]

        [grid.Buses[bus].V_magnitude_k = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ_k = JuMP.value.(δ[bus,!]) for bus in Nodes_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,0]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,0]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        [grid.Branches[branch].PowerFlow_ij_k = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,!]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji_k = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,!]) for branch in Branch_set]
        [grid.Branches[branch].losses_P_k = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij_k = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji_k = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q_k = 0 for branch in Branch_set]

        update_grid_tables!(grid)
    end

    return solution_status
end

function _solve_AC_SCOPF_TC!(grid ::PowerGrid,Redispatch=false,reference_gen=nothing,a=[],k=[])
    
end

function _solve_DC_SCOPF_TGC!(grid ::PowerGrid)
    
end


function _solve_AC_SCOPF_TGC!(grid ::PowerGrid)
    
end


include("SCOPF.jl")

@with_kw mutable struct SCOPF_Simulation_Settings
    # const available_preventive_actions = [:TransmssionSwitching,:BusbarSplitting,:SubstationsReconfiguration]
    # const available_corrective_actions = [:Redispatch, :TransmssionSwitching, :BusbarSplitting,:SubstationsReconfiguration]
    # const available_contingency_types = [:TC, :TGC] # transmission contingency, and transmission + generation contintgency
    preventive_actions = [:BusbarSplitting,:SubstationsReconfiguration]
    corrective_actions = [:TransmssionSwitching, :BusbarSplitting]
    formulation = :DCOPF # only :DCOPF formulation is available now
    method = :CCG # Column and Constraint Generation decomposition or :FullModel 
    contingency_type = :TC # can be :TransGen for considering transmission and generation contingencies
end

const default_simulation_settings = SCOPF_Simulation_Settings()

# TC: Transmission Contingencies -> Redispatch is optional
# TGC: Transmission and Generation Contingencies -> Redispatch is mandatory

function Generate_Contingencies(grid :: PowerGrid,  type ::Symbol)

    N_Line = length([grid.Branches[branch_id] for branch_id in keys(grid.Branches) if grid.Branches[branch_id].BranchType==0])
    N_Gen = length([grid.Generators[g] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual])
    N_Nodes = grid.N_bus;
    
    if type == :TC # usually used with Preventive-SCOPF
        k = N_Line+1
        a = ones(N_Nodes,N_Nodes,k)
        c = 1
        for branch_id in keys(grid.Branches)
            if grid.Branches[branch_id].BranchType == 0
                # doing it only for the physical non-auxilliary branches
                c += 1
                Fr_bus_ID = grid.Branches[branch_id].Fr_bus_ID
                To_bus_ID = grid.Branches[branch_id].To_bus_ID
                a[Fr_bus_ID,To_bus_ID,c] = 0
                a[To_bus_ID,Fr_bus_ID,c] = 0
            end
        end

        return a, k

    elseif type == :TGC # usually used with Corrective-SCOPF

        k = N_Line + N_Gen + 1;
        a_l = ones(N_Nodes,N_Nodes,k)
        a_g = ones(N_Gen,1,k)

        c = 1
        for branch_id in keys(grid.Branches)
            if grid.Branches[branch_id].BranchType == 0
                c += 1
                Fr_bus_ID = grid.Branches[branch_id].Fr_bus_ID
                To_bus_ID = grid.Branches[branch_id].To_bus_ID
                a_l[Fr_bus_ID,To_bus_ID,c] = 0
                a_l[To_bus_ID,Fr_bus_ID,c] = 0
            end
        end

        i = 1;
        for c in N_Line+2:k
            a_g[i,1,c] = 0;
            i = i+1;
        end

        return a_g,a_l,k
    end
end

function _solve_simulation!(grid ::PowerGrid,simulation_settings ::SCOPF_Simulation_Settings = default_simulation_settings)

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

    if simulation_settings.contingency_type == :Transmission
        a_l,K_ = Generate_Contingencies(System,:Transmission)
    elseif simulation_settings.contingency_type == :TransGen
        a_g,a_l,K_
    else
        println("Wrong contingency_type!! Aborting!!")
        return -1
    end

    contingecy_set = 1:K_
    # 1. Variables
    JuMP.@variable(m, grid.Buses[i,k].δ_min ≤ δ[i in Nodes_set, k in contingecy_set] ≤ grid.Buses[i,k].δ_max)

    if :Redispatch in simulation_settings.corrective_actions
        JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set, k in contingecy_set] ≤ grid.Generators[g].Pg_max)
    else
        JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    end

    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, k = contingecy_set; Set([i,j]) in Branch_nodes])

    # 2. Constraints

    # reference angle
    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[k = contingecy_set], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID, k] ==  0.0)
    end

    # Generator technical limits
    # if simulation_settings.contingency_type == :Transmission
    #     JuMP.@constraint(m,technical_limits[g in Gen_set],)
    # elseif simulation_settings.contingency_type == :TransGen

    # end
    
end

function _model_generator!(model,grid)
end


function solve_SCOBS_TC!(grid ::PowerGrid; Redispatch=false,reconf=[:pre,:post],splitting=[:pre, :post],reference_gen=nothing)
    if reconf == [:post] && splitting == [:pre]
        println("Wrong problem setting!")
        return -1
    elseif reconf == [:pre,:post] && length(splitting) == 1
        println("Wrong problem setting!")
        return -1
    elseif splitting == [:pre,:post] && reconf == [:post]
        println("Wrong problem setting!")
        return -1
    end

    a,c = Generate_Contingencies(grid,:TC)


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
    k_set = collect(1:c)

    # 1. Variables
    JuMP.@variable(m, grid.Buses[i].δ_min ≤ δ[i in Nodes_set,k in k_set] ≤ grid.Buses[i].δ_max)
    JuMP.@variable(m, pij[i = Nodes_set, j = Nodes_set, k in k_set; Set([i,j]) in Branch_nodes])
    
    if splitting == [:pre]
        JuMP.@variable(m, z_c[coupler = Coupler_set], Bin)
    else
        JuMP.@variable(m, z_c[coupler = Coupler_set, k in k_set], Bin)
    end

    if reconf == [:pre]
        JuMP.@variable(m, z_l[reconf = Reconf_set], Bin)
    else
        JuMP.@variable(m, z_l[reconf = Reconf_set, k in k_set], Bin)
    end

    if Redispatch
        JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set, k in k_set] ≤ grid.Generators[g].Pg_max)
        JuMP.@constraint(m,ramp_up[g in Gen_set,k in k_set[2:c]], p[g,k]-p[g,1] ≤ grid.Generators[g].Δ_up)
        JuMP.@constraint(m,ramp_down[g in Gen_set,k in k_set[2:c]], p[g,1]-p[g,k] ≤ grid.Generators[g].Δ_down)
    else
        JuMP.@variable(m, grid.Generators[g].Pg_min ≤ p[g in Gen_set] ≤ grid.Generators[g].Pg_max)
    end

    if Virtual_Gen_set != []
        # 2.0 virtual generator settings
        if splitting == [:pre]
            JuMP.@variable(m, p_vsc[g in Virtual_Gen_set])
            JuMP.@constraint(m,VSC_cap1[g in Virtual_Gen_set],(1-z_c[virtual_generator_couplers[g]])*grid.Generators[g].Pg_min ≤ p_vsc[g])
            JuMP.@constraint(m,VSC_cap2[g in Virtual_Gen_set], p_vsc[g] ≤ (1-z_c[virtual_generator_couplers[g]])*grid.Generators[g].Pg_max)
            JuMP.@constraint(m,VSC[l in Coupler_set],p_vsc[grid.Buses[grid.Branches[l].Fr_bus_ID].ConnectedGensIDs[1]] == -p_vsc[grid.Buses[grid.Branches[l].To_bus_ID].ConnectedGensIDs[1]])
        else
            JuMP.@variable(m, p_vsc[g in Virtual_Gen_set,k in k_set])
            JuMP.@constraint(m,VSC_cap1[g in Virtual_Gen_set,k in k_set],(1-z_c[virtual_generator_couplers[g],k])*grid.Generators[g].Pg_min ≤ p_vsc[g,k])
            JuMP.@constraint(m,VSC_cap2[g in Virtual_Gen_set,k in k_set], p_vsc[g,k] ≤ (1-z_c[virtual_generator_couplers[g],k])*grid.Generators[g].Pg_max)
            JuMP.@constraint(m,VSC[l in Coupler_set,k in k_set],p_vsc[grid.Buses[grid.Branches[l].Fr_bus_ID].ConnectedGensIDs[1],k]== -p_vsc[grid.Buses[grid.Branches[l].To_bus_ID].ConnectedGensIDs[1],k])
        end
    end
    

    # 2. Constraints

    if reference_gen !== nothing
        JuMP.@constraint(m, ReferenceAngle[k in k_set], δ[grid.Generators[Gen_set[reference_gen]].GenBus_ID,k] ==  0.0)
    end

    # Nodal balance constraint
    if Redispatch

        if Virtual_Gen_set != []

            if splitting == [:pre]
                # redispatch and only pre-contingency hybrid splitting
                JuMP.@constraint(m, Nodal_balance[i in Nodes_set, k in k_set],
                    sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p_vsc[g] for g in Virtual_Gen_set if grid.Generators[g].GenBus_ID == i) + sum(p[g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
            else
                # redispatch and pre+post -contingency hybrid splitting
                JuMP.@constraint(m, Nodal_balance[i in Nodes_set, k in k_set],
                    sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p_vsc[g,k] for g in Virtual_Gen_set if grid.Generators[g].GenBus_ID == i) + sum(p[g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
            end

        else
            # redispatch without virtual generator i.e. hybrid splitting
            JuMP.@constraint(m,Nodal_balance[i in Nodes_set,k in k_set],
                sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

        end
    else
        # no redispatch
        if Virtual_Gen_set != []

            if split == :pre
                # only pre-contingency hybrid splitting
                JuMP.@constraint(m, Nodal_balance[i in Nodes_set, k in k_set],
                    sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p_vsc[g] for g in Virtual_Gen_set if grid.Generators[g].GenBus_ID == i) + sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
            else
                # pre+post -contingency hybrid splitting
                JuMP.@constraint(m, Nodal_balance[i in Nodes_set, k in k_set],
                    sum(pij[i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(p_vsc[g,k] for g in Virtual_Gen_set if grid.Generators[g].GenBus_ID == i) + sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))
            end

        else
            # no virtual generator i.e. hybrid splitting
            JuMP.@constraint(m,Nodal_balance[i in Nodes_set,k in k_set],
                sum(pij[i,j,k] for j in Nodes_set if Set([i,j]) in Branch_nodes) == sum(p[g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) - sum(grid.Loads[d].Pd for d in keys(grid.Loads) if grid.Loads[d].LoadBus_ID == i))

        end

    end

    # 2.1 Active power flow accross transmission lines -> DCOPF equations
    JuMP.@constraint(m, pl[i in Nodes_set,j in Nodes_set,k in k_set; Set([i,j]) in Transmission_nodes],
        pij[i,j,k] == Sbase*(B[i,j])*(δ[i,k]-δ[j,k]))

    # 2.2 Opposite flow consistency
    JuMP.@constraint(m, pl_consistency[i in Nodes_set,j in Nodes_set,k in k_set; Set([i,j]) in Branch_nodes], 
        pij[i,j,k] == -pij[j,i,k])

    # 2.4 Thermal limits of transmission lines bounded by contingencies
    JuMP.@constraint(m, pl_rate[i in Nodes_set,j in Nodes_set,k in k_set; Set([i,j]) in Transmission_nodes],
        -a[i,j,k]*Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ pij[i,j,k] ≤ a[i,j,k]*Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)

    
    # Reconfiguration constraints
    if reconf == [:pre]
        # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
        JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set],sum(z_l[l] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
        
        # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes],
            δ[j,k]-(1-z_l[Reconf_dict[ Set([i,j])]]) ≤ δ[i,k]) 

        JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes],
            δ[i,k] ≤ δ[j,k] + (1-z_l[Reconf_dict[ Set([i,j])]]))

        # 2.7.1 Reconfiguration line capacity
        JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set, k in k_set; Set([i,j]) in Reconf_nodes], 
            -z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,k])
    
        JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set, k in k_set; Set([i,j]) in Reconf_nodes], 
            pij[i,j,k] ≤ z_l[Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)
    else
        # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
        JuMP.@constraint(m, no_circular_path_constraint[bus in aux_bus_set,k in k_set],sum(z_l[l,k] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
        
        # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes],
        δ[j,k]-(1-z_l[Reconf_dict[Set([i,j])],k]) ≤ δ[i,k]) 

        JuMP.@constraint(m, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes],
            δ[i,k] ≤ δ[j,k] + (1-z_l[Reconf_dict[ Set([i,j])],k]))
        
        # 2.7.1 Reconfiguration line capacity
        JuMP.@constraint(m,reconf_cap_1[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes], 
            -z_l[Reconf_dict[Set([i,j])],k] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ pij[i,j,k])
    
        JuMP.@constraint(m,reconf_cap_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Reconf_nodes], 
            pij[i,j,k] ≤ z_l[Reconf_dict[Set([i,j])],k] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)
    
    end
    
    # Splitting Constraints
    if splitting == [:pre]
        # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in k_set; Set([i,j]) in Coupler_nodes],
            δ[j,k] - (1-z_c[Coupler_dict[ Set([i,j])]]) ≤ δ[i,k] ) 

        JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes],
            δ[i,k] ≤ δ[j,k] + (1-z_c[Coupler_dict[ Set([i,j])]]) ) 

        # 2.7.2 Coupler line capacity
        JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes], 
            -z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,k])

        JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes], 
            pij[i,j,k] ≤ z_c[Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)
    else
        # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
        JuMP.@constraint(m, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in k_set; Set([i,j]) in Coupler_nodes],
            δ[j,k] - (1-z_c[Coupler_dict[ Set([i,j])],k]) ≤ δ[i,k] ) 

        JuMP.@constraint(m, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes],
            δ[i,k] ≤ δ[j,k] + (1-z_c[Coupler_dict[ Set([i,j])],k]) ) 

        # 2.7.2 Coupler line capacity
        JuMP.@constraint(m,coupler_cap_1[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes], 
            -z_c[Coupler_dict[Set([i,j])],k] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ pij[i,j,k])

        JuMP.@constraint(m,coupler_cap_2[i in Nodes_set, j in Nodes_set,k in k_set; Set([i,j]) in Coupler_nodes], 
            pij[i,j,k] ≤ z_c[Coupler_dict[Set([i,j])],k] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)
    end

    if splitting == [:post]
        JuMP.@constraint(m,pre_contingency_merge[coupler in Coupler_set], z_c[coupler, 1] == 1)
    end

    if reconf == [:post]
        JuMP.@constraint(m,pre_contingency_conf[l in Reconf_set], z_l[l, 1] == grid.Branches[l].GeneralSwitch.SwitchingStatus)
    end

    # 3. Objective
    if Redispatch
        JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g,1]+grid.Generators[g].C0 for g in Gen_set))
    else
        JuMP.@objective(m,Min,sum(grid.Generators[g].C1*p[g]+grid.Generators[g].C0 for g in Gen_set))
    end

    JuMP.optimize!(m)
    solution_status = raw_status(m)
    
    if solution_status != "Model was proven to be infeasible."

        # set new grid state
        grid.Operating_Cost = JuMP.objective_value(m)
        if Redispatch
            [grid.Generators[g].Pg = JuMP.value.(p[g,1]) for g in Gen_set]
            [grid.Generators[g].Qg = 0 for g in Gen_set]
            [grid.Generators[g].Pg_k = JuMP.value.(p[g,:]) for g in Gen_set]
            [grid.Generators[g].Qg_k = 0 for g in Gen_set]
        else
            [grid.Generators[g].Pg = JuMP.value.(p[g]) for g in Gen_set]
            [grid.Generators[g].Qg = 0 for g in Gen_set]
            [grid.Generators[g].Pg_k = JuMP.value.(p[g]) for g in Gen_set]
            [grid.Generators[g].Qg_k = 0 for g in Gen_set]
        end

        if Virtual_Gen_set != []
            if splitting == [:pre]
                [grid.Generators[g].Pg = JuMP.value.(p_vsc[g]) for g in Virtual_Gen_set]
                [grid.Generators[g].Qg = 0 for g in Virtual_Gen_set]
            else
                [grid.Generators[g].Pg = JuMP.value.(p_vsc[g,1]) for g in Virtual_Gen_set]
                [grid.Generators[g].Qg = 0 for g in Virtual_Gen_set]

                [grid.Generators[g].Pg_k = JuMP.value.(p_vsc[g,k]) for g in Virtual_Gen_set, k in k_set]
                [grid.Generators[g].Qg_k = 0 for g in Virtual_Gen_set]
            end
        end

        [grid.Buses[bus].V_magnitude = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ = JuMP.value.(δ[bus,1]) for bus in Nodes_set]

        [grid.Buses[bus].V_magnitude_k = 1 for bus in Nodes_set]
        [grid.Buses[bus].δ_k = JuMP.value.(δ[bus,k]) for bus in Nodes_set , k in k_set]

        [grid.Branches[branch].PowerFlow_ij = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,1]) for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,1]) for branch in Branch_set]
        [grid.Branches[branch].losses_P = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q = 0 for branch in Branch_set]

        [grid.Branches[branch].PowerFlow_ij_k = JuMP.value.(pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,k]) for branch in Branch_set, k in k_set]
        [grid.Branches[branch].PowerFlow_ji_k = JuMP.value.(pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,k]) for branch in Branch_set, k in k_set]
        [grid.Branches[branch].losses_P_k = grid.Branches[branch].PowerFlow_ij+grid.Branches[branch].PowerFlow_ji for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij_k = 0 for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji_k = 0 for branch in Branch_set]
        [grid.Branches[branch].losses_Q_k = 0 for branch in Branch_set]

        if reconf == [:pre]
            [grid.Branches[Reconf_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_l[Reconf_id]) for Reconf_id in Reconf_set]
        else
            [grid.Branches[Reconf_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_l[Reconf_id,1]) for Reconf_id in Reconf_set]
        end

        if splitting == [:pre]
            [grid.Branches[coupler_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_c[coupler_id]) for coupler_id in Coupler_set]
        else
            [grid.Branches[coupler_id].GeneralSwitch.SwitchingStatus = JuMP.value.(z_c[coupler_id,1]) for coupler_id in Coupler_set]
        end

        grid.Z_lines = JuMP.value.(z_l)
        grid.Z_coupler = JuMP.value.(z_c)

        for bus in aux_bus_set
            if reconf == [:pre]
                if sum(JuMP.value.(z_l[l]) for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) < 1
                    grid.Buses[bus].GeneralSwitch.SwitchingStatus = 0
                end
            else
                if sum(JuMP.value.(z_l[l,1]) for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) < 1
                    grid.Buses[bus].GeneralSwitch.SwitchingStatus = 0
                end
            end
        end

        reduce_grid!(grid)

        update_grid_tables!(grid)
    end

    return solution_status
end

function solve_SCOBS_TGC!(grid ::PowerGrid,reconf=[:Pre,:Post],splitting=[:pre, :Post])
    
end