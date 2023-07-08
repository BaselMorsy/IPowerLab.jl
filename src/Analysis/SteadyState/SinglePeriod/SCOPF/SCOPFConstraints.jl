using JuMP
using Parameters

@with_kw mutable struct SCOPF_SimulationSettings
    ac_grid_model = :DCOPF
    dc_grid_model = :Linear
    converter_model = :Linear_lossless
    dynamic_converter_control = true
    converter_modularization = :discrete

    transmission_switching = [] # [:pre] , [:post] , [:pre,:post]
    substation_switching = Dict("splitting" => [:post], "reconf" => [:pre])
    max_transmission_switching = Inf
    max_substation_reconf = Inf
    max_busbar_splitting = Inf
    
    redispatch = false

    contingency_type = :Transmission_Contingencies # :TransGen_Contingencies
    
    NLP_solver = Ipopt.Optimizer
    MILP_solver = Gurobi.Optimizer
    Meta_solver = nothing # => can be :CCG , :RL , :RL_OPT
end

@with_kw mutable struct SCOPF_Prerequisites
    nodes_set
    aux_bus_set
    branch_nodes
    gen_set
    load_set
    ac_virtual_gen_set
    dc_virtual_gen_set
    all_gen_set
    B_matrix
    base_MVA
    unswitched_Transmission_nodes
    branch_dictionary
    switched_Transmission_nodes
    reconf_nodes
    reconf_dict
    coupler_nodes
    coupler_dict
    default_off_reconf
    coupler_set
    converter_set
    dc_Nodes_set
    dc_Transmission_nodes
    dc_branch_dictionary
    b2b_gen_set
    b2b_coupler_dict
    switched_transmission_set
    reconf_set
    a_g
    a_l
    a_l_dc_branch
    a_l_dc_link
    k
    reference_node = nothing
end

function single_period_SCOPF_variable_initialization!(model::Model, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites)
    """
    Variable initialization for SCOPF simulation, based on simulation settings `simulation_settings` and prerequisites `prerequisites_data`
    """
    # variable initialization
    if simulation_settings.ac_grid_model == :DCOPF
        JuMP.@variable(model, δ[i in prerequisites_data.nodes_set, k in prerequisites_data.k])
        JuMP.@variable(model, pij[i in prerequisites_data.nodes_set, j in prerequisites_data.nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.branch_nodes])
        JuMP.@variable(model, p[g in prerequisites_data.gen_set, k in prerequisites_data.k])
    elseif simulation_settings.ac_grid_model == :ACOPF
        JuMP.@variable(model, v[i in prerequisites_data.nodes_set, k in prerequisites_data.k])
        JuMP.@variable(model, δ[i in prerequisites_data.nodes_set, k in prerequisites_data.k])
        JuMP.@variable(model, pij[i in prerequisites_data.nodes_set, j in prerequisites_data.nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.branch_nodes])
        JuMP.@variable(model, qij[i in prerequisites_data.nodes_set, j in prerequisites_data.nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.branch_nodes])

        JuMP.@variable(model, p[g in prerequisites_data.gen_set, k in prerequisites_data.k])
        JuMP.@variable(model, q[g in prerequisites_data.gen_set, k in prerequisites_data.k])
    else
        error("Unidentified grid model setting ($simulation_settings.ac_grid_model)")
        return -1
    end

    if simulation_settings.dc_grid_model == :Linear
        JuMP.@variable(model, pij_dc[i in prerequisites_data.dc_Nodes_set, j in prerequisites_data.dc_Nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.dc_Transmission_nodes])
    elseif simulation_settings.dc_grid_model == :NonLinear
        JuMP.@variable(model, v_dc[i in prerequisites_data.dc_Nodes_set, k in prerequisites_data.k])
        JuMP.@variable(model, pij_dc[i in prerequisites_data.dc_Nodes_set, j in prerequisites_data.dc_Nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.dc_Transmission_nodes])
    else
        error("Unidentified grid model setting ($simulation_settings.dc_grid_model)")
        return -1
    end

    # converter variables
    JuMP.@variable(model, p_conv[g in collect(union(Set(prerequisites_data.ac_virtual_gen_set), Set(prerequisites_data.dc_virtual_gen_set))), k in prerequisites_data.k])


    # Switching variables
    if prerequisites_data.switched_transmission_set != []
        JuMP.@variable(model, z[l in prerequisites_data.switched_transmission_set, k in prerequisites_data.k], Bin)
    end

    if simulation_settings.substation_switching["reconf"] == [:post] && simulation_settings.substation_switching["splitting"] == [:pre]
        error("Wrong problem setting!")
        return -1
    elseif simulation_settings.substation_switching["reconf"] == [:pre,:post] && length(simulation_settings.substation_switching["splitting"]) == 1
        error("Wrong problem setting!")
        return -1
    elseif simulation_settings.substation_switching["splitting"] == [:pre,:post] && simulation_settings.substation_switching["reconf"] == [:post]
        error("Wrong problem setting!")
        return -1
    end

    if prerequisites_data.reconf_set != []
        JuMP.@variable(model, z_l[l in prerequisites_data.reconf_set, k in prerequisites_data.k], Bin)
    end

    if prerequisites_data.coupler_set != []
        JuMP.@variable(model, z_c[c in prerequisites_data.coupler_set, k in prerequisites_data.k], Bin)
    end
end

function single_period_SCOPF_LoadShedding_variable_initialization!(model::Model, prerequisites_data::SCOPF_Prerequisites)
    JuMP.@variable(model, p_ls[d in prerequisites_data.load_set, k in prerequisites_data.k])
end

function single_period_SC_nodal_balance_ac_node!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Nodes_set = prereq.nodes_set
    Branch_nodes = prereq.branch_nodes
    load_set = prereq.load_set
    Gen_set = prereq.gen_set
    ac_virtual_gen_set = prereq.ac_virtual_gen_set

    if formulation == :ACOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set, k in prereq.k],
            sum(model[:pij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_conv][g,k] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
        
        # REACTIVE NODAL BALANCE
        JuMP.@constraint(model, Qnodal[i in Nodes_set,k in prereq.k],
            sum(model[:qij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:q][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i)
                 - sum(grid.Loads[l].Qd for l in load_set if grid.Loads[l].LoadBus_ID == i))
    elseif formulation == :DCOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set,k in prereq.k],
            sum(model[:pij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_conv][g,k] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_nodal_balance_LS_ac_node!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Nodes_set = prereq.nodes_set
    Branch_nodes = prereq.branch_nodes
    load_set = prereq.load_set
    Gen_set = prereq.gen_set
    ac_virtual_gen_set = prereq.ac_virtual_gen_set

    if formulation == :ACOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set, k in prereq.k],
            sum(model[:pij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_ls][l,k] for l in load_set if grid.Loads[l].LoadBus_ID == i)
                 + sum(model[:p_conv][g,k] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
        
        # REACTIVE NODAL BALANCE
        JuMP.@constraint(model, Qnodal[i in Nodes_set,k in prereq.k],
            sum(model[:qij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:q][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i)
                 - sum(grid.Loads[l].Qd for l in load_set if grid.Loads[l].LoadBus_ID == i))
    elseif formulation == :DCOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set,k in prereq.k],
            sum(model[:pij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_conv][g,k] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i)
                + sum(model[:p_ls][l,k] for l in load_set if grid.Loads[l].LoadBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_powerflow_ac_branch!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Nodes_set = prereq.nodes_set
    Branch_nodes = prereq.branch_nodes
    unswitched_Transmission_nodes = prereq.unswitched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    Sbase = prereq.base_MVA
    B = prereq.B_matrix

    if formulation == :ACOPF
        # ACTIVE POWER THROUGH LINE i-j
        JuMP.@NLconstraint(model, p_line[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i,j]) in unswitched_Transmission_nodes],
            model[:pij][i,j,k] ==  prereq.a_l[i,j,k]*Sbase*(model[:v][i,k]*model[:v][j,k]*(G[i,j]*cos(model[:δ][i,k]-model[:δ][j,k])+B[i,j]*sin(model[:δ][i,k]-model[:δ][j,k])) -(model[:v][i,k]^2)*G[i,j]))

        # REACTIVE POWER THROUGH LINE i-j
        JuMP.@NLconstraint(model, q_line[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i,j]) in unswitched_Transmission_nodes],
            model[:pij][i,j,k] ==  prereq.a_l[i,j,k]*Sbase*(model[:v][i,k]*model[:v][j,k]*(G[i,j]*sin(model[:δ][i,k]-model[:δ][j,k]) - B[i,j]*cos(model[:δ][i,k]-model[:δ][j,k])) +(model[:v][i,k]^2)*(B[i,j] - grid.Branches[branch_dictionary[Set([i,j])]].b/2)))

    elseif formulation == :DCOPF
        # ACTIVE POWER THROUGH LINE i-j
        JuMP.@constraint(model,pl[i in Nodes_set,j in Nodes_set, k in prereq.k; Set([i,j]) in unswitched_Transmission_nodes],
            model[:pij][i,j,k] == Sbase*(B[i,j])*(model[:δ][i,k]-model[:δ][j,k])*prereq.a_l[i,j,k])
        
        # OPPOSITE FLOW CONSISTENCY
        JuMP.@constraint(model, pl_consistency[i in Nodes_set,j in Nodes_set, k in prereq.k; Set([i,j]) in Branch_nodes], 
        model[:pij][i,j,k] == -model[:pij][j,i,k])
    else
        error("Requested formulation ($formulation) is not implemented.")
    end  
end

function single_period_SC_transmission_capacity_limits_ac_branch!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)

    Nodes_set = prereq.nodes_set
    Sbase = prereq.base_MVA
    branch_dictionary = prereq.branch_dictionary
    unswitched_Transmission_nodes = prereq.unswitched_Transmission_nodes

    if formulation == :ACOPF
        # LINE CAPACITY
        JuMP.@NLconstraint(model, Smax[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i,j]) in unswitched_Transmission_nodes],
            (model[:pij][i,j,k])^2 + (model[:qij][i,j,k])^2 ≤ prereq.a_l[i,j,k] * ((Sbase)^2)*(grid.Branches[branch_dictionary[Set([i,j])]].rating)^2)
    elseif formulation == :DCOPF
        # LINE CAPACITY
        JuMP.@constraint(model,pl_rate[i in Nodes_set,j in Nodes_set, k in prereq.k; Set([i,j]) in unswitched_Transmission_nodes],
            -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*prereq.a_l[i,j,k] ≤ model[:pij][i,j,k] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*prereq.a_l[i,j,k])
    else
        error("Requested formulation ($formulation) is not implemented.")
    end  
    
end

function single_period_SC_switched_transmission_capacity_limits_ac_branch!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model)

    Sbase = prereq.base_MVA
    Nodes_set = prereq.nodes_set
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    # LINE CAPACITY
    JuMP.@constraint(model, pl_rate_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        -Sbase * grid.Branches[branch_dictionary[Set([i, j])]].rating * model[:z][branch_dictionary[Set([i, j])],k] * prereq.a_l[i,j,k] ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, pl_rate_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        model[:pij][i, j, k] ≤ Sbase * grid.Branches[branch_dictionary[Set([i, j])]].rating * model[:z][branch_dictionary[Set([i, j])],k] * prereq.a_l[i,j,k])

end

function single_period_SC_SG_switched_transmission_capacity_limits_ac_branch!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model,z)
    # SG: static grid
    Sbase = prereq.base_MVA
    Nodes_set = prereq.nodes_set
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    # LINE CAPACITY
    JuMP.@constraint(model,pl_rate_1[i in Nodes_set,j in Nodes_set, k in prereq.k; Set([i,j]) in switched_Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])],k] * prereq.a_l[i,j,k] ≤ model[:pij][i,j,k] )

    JuMP.@constraint(model,pl_rate_2[i in Nodes_set,j in Nodes_set, k in prereq.k; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j,k] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])],k] * prereq.a_l[i,j,k])

end

function single_period_SC_switched_powerflow_ac_branch!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model; max_op=Inf)

    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    Nodes_set = prereq.nodes_set
    branch_dictionary = prereq.branch_dictionary
    B = prereq.B_matrix
    Sbase = prereq.base_MVA
    M = grid.S_base * 100
    JuMP.@constraint(model, pl_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i,k] - model[:δ][j,k]) * prereq.a_l[i,j,k] ≤ (1 - model[:z][branch_dictionary[Set([i, j])],k]) * M)

    JuMP.@constraint(model, pl_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i,k] - model[:δ][j,k]) * prereq.a_l[i,j,k] ≥ -(1 - model[:z][branch_dictionary[Set([i, j])],k]) * M)

    if max_op !== Inf
        JuMP.@constraint(model, maximum_allowed_switching_actions[k in prereq.k], sum(1 - model[:z][branch_dictionary[l],k] for l in switched_Transmission_nodes) ≤ max_op)
    end
end

function single_period_SC_SG_switched_powerflow_ac_branch!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model, z)

    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    Nodes_set = prereq.nodes_set
    branch_dictionary = prereq.branch_dictionary
    B = prereq.B_matrix
    Sbase = prereq.base_MVA
    M = Sbase * 100
    JuMP.@constraint(model, pl_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i, k] - model[:δ][j, k]) * prereq.a_l[i,j,k] ≤ (1 - z[branch_dictionary[Set([i, j])], k]) * M)

    JuMP.@constraint(model, pl_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in switched_Transmission_nodes],
        model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i, k] - model[:δ][j, k]) * prereq.a_l[i,j,k] ≥ -(1 - z[branch_dictionary[Set([i, j])], k]) * M)
end

function single_period_SC_reconf_split_constraints_ac_grid!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model; max_reconf=Inf, max_splitting=Inf)

    aux_bus_set = prereq.aux_bus_set
    Reconf_dict = prereq.reconf_dict
    Nodes_set = prereq.nodes_set
    Reconf_nodes = prereq.reconf_nodes
    Coupler_nodes = prereq.coupler_nodes
    Coupler_dict = prereq.coupler_dict
    Sbase = prereq.base_MVA
    default_off_reconf = prereq.default_off_reconf
    Coupler_set = prereq.coupler_set

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(model, no_circular_path_constraint[bus in aux_bus_set, k in prereq.k], sum(model[:z_l][l, k] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
    M_δ = 2 * π

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][j, k] - M_δ * (1 - model[:z_l][Reconf_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - model[:z_l][Reconf_dict[Set([i, j])], k]))

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][j, k] - M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]))

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(model, reconf_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        -model[:z_l][Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, reconf_cap_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:pij][i, j, k] ≤ model[:z_l][Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating)

    # 2.7.2 coupler capacity
    JuMP.@constraint(model, coupler_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        -model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, coupler_cap_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:pij][i, j, k] ≤ model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating)

    if max_reconf !== Inf
        JuMP.@constraint(model, max_reconf_actions[k in prereq.k], sum(model[:z_l][r,k] for r in default_off_reconf) ≤ max_reconf)
    end

    if max_splitting !== Inf
        JuMP.@constraint(model, max_splitting_actions[k in prereq.k], sum(model[:z_c][c,k] for c in Coupler_set) ≤ max_splitting)
    end

end

function single_period_SC_reconf_split_SG_constraints_ac_grid!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model, z_l, z_c)

    aux_bus_set = prereq.aux_bus_set
    Reconf_dict = prereq.reconf_dict
    Nodes_set = prereq.nodes_set
    Reconf_nodes = prereq.reconf_nodes
    Coupler_nodes = prereq.coupler_nodes
    Coupler_dict = prereq.coupler_dict
    Sbase = prereq.base_MVA
    default_off_reconf = prereq.default_off_reconf
    Coupler_set = prereq.coupler_set

    M_δ = 2 * π

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][j, k] - M_δ * (1 - z_l[Reconf_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - z_l[Reconf_dict[Set([i, j])], k]))

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][j, k] - M_δ * (1 - z_c[Coupler_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - z_c[Coupler_dict[Set([i, j])], k]))

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(model, reconf_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        -z_l[Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        model[:pij][i, j] ≤ z_l[Reconf_dict[Set([i, j])]] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(model, coupler_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        -z_c[Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, coupler_cap_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:pij][i, j, k] ≤ z_c[Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating)

end

function single_period_SC_reconf_SG_constraints_ac_grid!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model, z_l)

    aux_bus_set = prereq.aux_bus_set
    Reconf_dict = prereq.reconf_dict
    Nodes_set = prereq.nodes_set
    Reconf_nodes = prereq.reconf_nodes
    Coupler_nodes = prereq.coupler_nodes
    Coupler_dict = prereq.coupler_dict
    Sbase = prereq.base_MVA
    default_off_reconf = prereq.default_off_reconf
    Coupler_set = prereq.coupler_set

    M_δ = 2 * π

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][j, k] - M_δ * (1 - z_l[Reconf_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - z_l[Reconf_dict[Set([i, j])], k]))

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][j, k] - M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

    JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]))

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(model, reconf_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Reconf_nodes],
        -z_l[Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        model[:pij][i, j] ≤ z_l[Reconf_dict[Set([i, j])]] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(model, coupler_cap_1[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        -model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

    JuMP.@constraint(model, coupler_cap_2[i in Nodes_set, j in Nodes_set, k in prereq.k; Set([i, j]) in Coupler_nodes],
        model[:pij][i, j, k] ≤ model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating)

end

function single_period_SC_generator_limits_ac_grid!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Gen_set = prereq.gen_set

    if formulation == :ACOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Pg_min * prereq.a_g[g,k] ≤ model[:p][g,k] ≤ grid.Generators[g].Pg_max * prereq.a_g[g,k])
        JuMP.@constraint(model, gen_reactive_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Qg_min * prereq.a_g[g,k] ≤ model[:q][g,k] ≤ grid.Generators[g].Qg_max * prereq.a_g[g,k])
    elseif formulation == :DCOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Pg_min * prereq.a_g[g,k] ≤ model[:p][g,k] ≤ grid.Generators[g].Pg_max * prereq.a_g[g,k])
        JuMP.@constraint(model, total_balance[k in prereq.k], sum(model[:p][g,k] for g in Gen_set) == sum(grid.Loads[l].Pd for l in prereq.load_set))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_generator_limits_LS_ac_grid!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model, formulation::Symbol)

    Gen_set = prereq.gen_set

    if formulation == :ACOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Pg_min * prereq.a_g[g, k] ≤ model[:p][g, k] ≤ grid.Generators[g].Pg_max * prereq.a_g[g, k])
        JuMP.@constraint(model, gen_reactive_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Qg_min * prereq.a_g[g, k] ≤ model[:q][g, k] ≤ grid.Generators[g].Qg_max * prereq.a_g[g, k])
    elseif formulation == :DCOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set, k in prereq.k], grid.Generators[g].Pg_min * prereq.a_g[g, k] ≤ model[:p][g, k] ≤ grid.Generators[g].Pg_max * prereq.a_g[g, k])
        JuMP.@constraint(model, total_balance[k in prereq.k], sum(model[:p][g, k] for g in Gen_set) == sum(grid.Loads[l].Pd for l in prereq.load_set) - sum(model[:p_ls][l,k] for l in load_set if grid.Loads[l].LoadBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_voltage_limits_ac_grid!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    if formulation == :ACOPF
        Nodes_set = prereq.nodes_set
        JuMP.@constraint(model, voltage_limits[i in Nodes_set,k in prereq.k], grid.Buses[i].V_min ≤ model[:v][i,k] ≤ grid.Buses[i].V_max)
    end
end

function single_period_SC_angle_limits_ac_grid!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model)
    Nodes_set = prereq.nodes_set
    JuMP.@constraint(model, angle_limits[i in Nodes_set, k in prereq.k], grid.Buses[i].δ_min ≤ model[:δ][i,k] ≤ grid.Buses[i].δ_max)
    if !isnothing(prereq.reference_node)
        JuMP.@constraint(model, reference_node[k in prereq.k], model[:δ][prereq.reference_node,k] == 0)
    end
end

function single_period_SC_converter_flow!(prereq::SCOPF_Prerequisites, grid::PowerGrid, model::Model, formulation::Symbol; modularization = :discrete)
    converter_set = prereq.converter_set
    b2b_gen_set = prereq.b2b_gen_set
    b2b_coupler_dict = prereq.b2b_coupler_dict
    ac_virtual_gen_set = prereq.ac_virtual_gen_set
    dc_virtual_gen_set = prereq.dc_virtual_gen_set

    if formulation == :Nonlinear
        error("NonLinear formulation for converters is not implemented yet")
    elseif formulation == :Linear_lossless
        if grid.N_conv_duplets == 0
            modularization = :discrete
        end

        if modularization == :discrete
            JuMP.@constraint(model, conv_power_flow[c in converter_set, k in prereq.k], model[:p_conv][grid.Converters[c].gen_ac_id, k] + model[:p_conv][grid.Converters[c].gen_dc_id, k] == 0)
        elseif modularization == :continuous
            duplets = grid.Converter_Duplets
            convs = grid.Converters
            JuMP.@constraint(model, continous_rating_ac_side[d in collect(keys(duplets)), k in prereq.k], grid.Converters[duplets[d][1]].rate ≤ model[:p_conv][convs[duplets[d][1]].gen_ac_id,k] + model[:p_conv][convs[duplets[d][2]].gen_ac_id,k] ≤ grid.Converters[duplets[d][1]].rate)
            JuMP.@constraint(model, continous_rating_dc_side[d in collect(keys(duplets)), k in prereq.k], grid.Converters[duplets[d][1]].rate ≤ model[:p_conv][convs[duplets[d][1]].gen_dc_id,k] + model[:p_conv][convs[duplets[d][2]].gen_dc_id,k] ≤ grid.Converters[duplets[d][1]].rate)
            JuMP.@constraint(model, conv_power_flow_no_loss[c in converter_set, k in prereq.k], model[:p_conv][grid.Converters[c].gen_ac_id, k] + model[:p_conv][grid.Converters[c].gen_dc_id, k] == 0)
            JuMP.@constraint(model, conv_power_flow[d in collect(keys(duplets)), k in prereq.k], model[:p_conv][convs[duplets[d][1]].gen_ac_id,k] + model[:p_conv][convs[duplets[d][2]].gen_ac_id,k] + model[:p_conv][convs[duplets[d][1]].gen_dc_id,k] + model[:p_conv][convs[duplets[d][2]].gen_dc_id,k] == 0)
        else
            error("Undefined modularization technique $(modularization)")
        end
    elseif formulation == :Linear_lossy
        error("Linear_lossy formulation for converters is not implemented yet")
        # JuMP.@constraint(model, conv_power[c in converter_set], p_conv[grid.Converters[c].gen_ac_id] + p_conv[grid.Converters[c].gen_dc_id] == losses)
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_converter_constraints!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    converter_set = prereq.converter_set
    b2b_gen_set = prereq.b2b_gen_set
    b2b_coupler_dict = prereq.b2b_coupler_dict
    ac_virtual_gen_set = prereq.ac_virtual_gen_set
    dc_virtual_gen_set = prereq.dc_virtual_gen_set
    
    if formulation == :Nonlinear
        error("NonLinear formulation for converters is not implemented yet")
    elseif formulation == :Linear_lossless
        # JuMP.@constraint(model, conv_power[c in converter_set,k in prereq.k], model[:p_conv][grid.Converters[c].gen_ac_id,k] + model[:p_conv][grid.Converters[c].gen_dc_id,k] == 0)
        JuMP.@constraint(model, dc_link_power[link_id in keys(grid.DCLinks), k in prereq.k], model[:p_conv][grid.DCLinks[link_id].Fr_gen_ID,k] + model[:p_conv][grid.DCLinks[link_id].To_gen_ID,k] == 0)
        JuMP.@constraint(model, conv_cap_ac_side[g in ac_virtual_gen_set, k in prereq.k], grid.Generators[g].Pg_min ≤ model[:p_conv][g,k] ≤ grid.Generators[g].Pg_max)
        JuMP.@constraint(model, conv_cap_dc_side[g in dc_virtual_gen_set, k in prereq.k], grid.Generators[g].Pg_min ≤ model[:p_conv][g,k] ≤ grid.Generators[g].Pg_max)
        if b2b_gen_set !== []
            # SOFT SPLITTING CONSTRAINTS
            JuMP.@constraint(model,VSC_cap1[g in b2b_gen_set, k in prereq.k],(1-model[:z_c][b2b_coupler_dict[g],k])*grid.Generators[g].Pg_min ≤ model[:p_conv][g,k])
            JuMP.@constraint(model,VSC_cap2[g in b2b_gen_set, k in prereq.k], model[:p_conv][g,k] ≤ (1-model[:z_c][b2b_coupler_dict[g],k])*grid.Generators[g].Pg_max)
        end

    elseif formulation == :Linear_lossy
        error("Linear_lossy formulation for converters is not implemented yet")
        # JuMP.@constraint(model, conv_power[c in converter_set], p_conv[grid.Converters[c].gen_ac_id] + p_conv[grid.Converters[c].gen_dc_id] == losses)
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_nodal_balance_dc_node!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes
    dc_virtual_gen_set = prereq.dc_virtual_gen_set
    
    if formulation == :NonLinear
        error("NonLinear nodal balance is not implemented yet for DC grids")
    elseif formulation == :Linear
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal_dc[i in dc_Nodes_set,k in prereq.k],
            sum(model[:pij_dc][i,j,k] for j = dc_Nodes_set if Set([i,j]) in dc_Transmission_nodes) == sum(model[:p_conv][g,k] for g in dc_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_SC_transmission_capacity_limits_dc_grid!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Sbase = prereq.base_MVA
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes
    dc_branch_dictionary = prereq.dc_branch_dictionary

    if formulation == :Linear
        # LINE CAPACITY
        JuMP.@constraint(model,pl_rate_dc[i in dc_Nodes_set,j in dc_Nodes_set,k in prereq.k; Set([i,j]) in dc_Transmission_nodes],
            -Sbase*grid.DCBranches[dc_branch_dictionary[Set([i,j])]].rating ≤ model[:pij_dc][i,j,k] ≤ Sbase*grid.DCBranches[dc_branch_dictionary[Set([i,j])]].rating)
    elseif formulation == :NonLinear
        error("NonLinear flow limits in DC grid is not implemented yet")
    else
        error("Requested formulation ($formulation) is not implemented")
    end
end

function single_period_SC_powerflow_dc_branch!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes

    if formulation == :Linear
        JuMP.@constraint(model,dc_Pij[i in dc_Nodes_set, j in dc_Nodes_set,k in prereq.k;  Set([i,j]) in dc_Transmission_nodes], model[:pij_dc][i,j,k] + model[:pij_dc][j,i,k] == 0 )
    elseif formulation == :NonLinear
        error("NonLinear power flow formulation for DC grid is not implemented yet")
    else
        error("Requested formulation ($formulation) is not implemented")
    end
end

function single_period_SC_voltage_limits_dc_grid!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    if formulation == :NonLinear
        dc_Nodes_set = prereq.dc_Nodes_set
        JuMP.@constraint(model, dc_voltage[i in dc_Nodes_set,k in prereq.k], grid.DCBuses[i].V_min ≤ model[:v_dc][i,k] ≤ grid.DCBuses[i].V_max)
    end
end

function single_period_SC_objective!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, transmission_switching ::Bool, substation_switching ::Bool)
    
    if transmission_switching && substation_switching
        Coupler_set = prereq.coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z_c][c,k] for c in Coupler_set , k in prereq.k)*1e-4+sum(1-model[:z][l,k] for l in prereq.switched_transmission_set, k in prereq.k)*1e-4)
    elseif transmission_switching && ! substation_switching
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z][l,k] for l in prereq.switched_transmission_set, k in prereq.k)*1e-4)
    elseif substation_switching && ! transmission_switching
        Coupler_set = prereq.coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z_c][c,k] for c in Coupler_set, k in prereq.k)*1e-4)
    else
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set))
    end

end

function single_period_SC_LS_objective!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, transmission_switching ::Bool, substation_switching ::Bool;LS_cost=1000)
    
    if transmission_switching && substation_switching
        Coupler_set = prereq.coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z_c][c,k] for c in Coupler_set , k in prereq.k)*1e-4 + sum(1-model[:z][l,k] for l in prereq.switched_transmission_set, k in prereq.k)*1e-4) + sum(model[:p_ls][l,k] for c in prereq.load_set , k in prereq.k)*LS_cost
    elseif transmission_switching && ! substation_switching
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z][l,k] for l in prereq.switched_transmission_set, k in prereq.k)*1e-4) + sum(model[:p_ls][l,k] for c in prereq.load_set , k in prereq.k)*LS_cost
    elseif substation_switching && ! transmission_switching
        Coupler_set = prereq.coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)+sum(1-model[:z_c][c,k] for c in Coupler_set, k in prereq.k)*1e-4) + sum(model[:p_ls][l,k] for c in prereq.load_set , k in prereq.k)*LS_cost
    else
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g,1]+grid.Generators[g].C0 for g in prereq.gen_set)) + sum(model[:p_ls][l,k] for c in prereq.load_set , k in prereq.k)*LS_cost
    end

end

function single_period_SC_corrective_measures!(prereq ::SCOPF_Prerequisites, grid ::PowerGrid, model ::Model, simulation_settings ::SCOPF_SimulationSettings)
    
    non_base_contingencies = collect(setdiff(Set(prereq.k),Set([1])))
    # Corrective redispatch constraints
    if ! simulation_settings.redispatch
        JuMP.@constraint(model,no_redispatch[g in prereq.gen_set, k in non_base_contingencies], model[:p][g,1] == model[:p][g,k])
    else
        JuMP.@constraint(model,Up_redispatch[g in prereq.gen_set, k in non_base_contingencies], model[:p][g,1] - model[:p][g,k] ≤ grid.Generators[g].Δ_down)
        JuMP.@constraint(model,Down_redispatch[g in prereq.gen_set, k in non_base_contingencies], model[:p][g,k] - model[:p][g,1] ≤ grid.Generators[g].Δ_up)
    end

    # Dynamic converter control constraints
    if ! simulation_settings.dynamic_converter_control
        JuMP.@constraint(model,static_converters[g in collect(union(Set(prereq.ac_virtual_gen_set), Set(prereq.dc_virtual_gen_set))), k in non_base_contingencies], model[:p_conv][g,1] == model[:p_conv][g,k])
    end

    # Corrective transmission switching constraints
    if simulation_settings.transmission_switching == [:post]
        JuMP.@constraint(model, static_grid_pre_contingency[l in prereq.switched_transmission_set], model[:z][l,1] == 1)
    elseif simulation_settings.transmission_switching == [:pre]
        JuMP.@constraint(model, static_grid_post_contingency[l in prereq.switched_transmission_set,k in non_base_contingencies], model[:z][l,k] == model[:z][l,1])
    end

    # Corrective substation switching constraints
    if simulation_settings.substation_switching["reconf"] == [:pre]
        JuMP.@constraint(model, static_substation_post_contingency[l in prereq.reconf_set, k in non_base_contingencies], model[:z_l][l,k] == model[:z_l][l,1])
    end

    if simulation_settings.substation_switching["splitting"] == [:pre]
        JuMP.@constraint(model, static_coupler_post_contingency[l in prereq.coupler_set, k in non_base_contingencies], model[:z_c][l,k] == model[:z_c][l,1])
    elseif simulation_settings.substation_switching["splitting"] == [:post]
        JuMP.@constraint(model, static_coupler_pre_contingency[l in prereq.coupler_set], model[:z_c][l,1] == 1)
    elseif simulation_settings.substation_switching["splitting"] == []
        JuMP.@constraint(model, static_coupler[l in prereq.coupler_set, k in prereq.k], model[:z_c][l,k] == 1)
    end
end

