using JuMP
using Parameters

@with_kw mutable struct OPF_SimulationSettings
    ac_grid_model = :DCOPF
    transmission_switching = true
    substation_switching = true
    max_transmission_switching = Inf
    max_substation_reconf = Inf
    max_busbar_splitting = Inf
    converter_model = :Linear_lossless
    dc_grid_model = :Linear
    NLP_solver = Ipopt.Optimizer
    MILP_solver = Gurobi.Optimizer
end

@with_kw mutable struct OPF_Prerequisites
    Nodes_set
    aux_bus_set
    Branch_nodes
    Gen_set
    load_set
    ac_virtual_gen_set
    dc_virtual_gen_set
    all_gen_set
    B
    Sbase
    unswitched_Transmission_nodes
    branch_dictionary
    switched_Transmission_nodes
    Reconf_nodes
    Reconf_dict
    Coupler_nodes
    Coupler_dict
    default_off_reconf
    Coupler_set
    converter_set
    dc_Nodes_set
    dc_Transmission_nodes
    dc_branch_dictionary
    b2b_gen_set
    b2b_coupler_dict
    switched_transmission_set
    reconf_set
    reference_node = nothing
end

function single_period_OPF_variable_initialization!(model ::Model, simulation_settings ::OPF_SimulationSettings, prerequisites_data ::OPF_Prerequisites; ideal_nodes = [])
    # variable initialization
    if simulation_settings.ac_grid_model == :DCOPF
        JuMP.@variable(model,δ[i in prerequisites_data.Nodes_set])

        JuMP.@constraint(model, fixed_δ[i in ideal_nodes], δ[i] == 0)

        JuMP.@variable(model,pij[i in prerequisites_data.Nodes_set,j in prerequisites_data.Nodes_set; Set([i,j]) in prerequisites_data.Branch_nodes])
        JuMP.@variable(model,p[g in prerequisites_data.Gen_set])
    elseif simulation_settings.ac_grid_model == :ACOPF
        JuMP.@variable(model,v[i in prerequisites_data.Nodes_set])
        JuMP.@variable(model,δ[i in prerequisites_data.Nodes_set])

        JuMP.@constraint(model, fixed_δ[i in ideal_nodes], δ[i] == 0)
        JuMP.@constraint(model, fixed_v[i in ideal_nodes], v[i] == 1)

        JuMP.@variable(model,pij[i in prerequisites_data.Nodes_set,j in prerequisites_data.Nodes_set; Set([i,j]) in prerequisites_data.Branch_nodes])
        JuMP.@variable(model,qij[i in prerequisites_data.Nodes_set,j in prerequisites_data.Nodes_set; Set([i,j]) in prerequisites_data.Branch_nodes])
        JuMP.@variable(model,p[g in prerequisites_data.Gen_set])
        JuMP.@variable(model,q[g in prerequisites_data.Gen_set])
    else
        error("Unidentified grid model setting ($simulation_settings.ac_grid_model)")
        return -1
    end

    if simulation_settings.dc_grid_model == :Linear
        JuMP.@variable(model,pij_dc[i in prerequisites_data.dc_Nodes_set,j in prerequisites_data.dc_Nodes_set; Set([i,j]) in prerequisites_data.dc_Transmission_nodes])
    elseif simulation_settings.dc_grid_model == :NonLinear
        JuMP.@variable(model,v_dc[i in prerequisites_data.dc_Nodes_set])
        JuMP.@variable(model,pij_dc[i in prerequisites_data.dc_Nodes_set,j in prerequisites_data.dc_Nodes_set; Set([i,j]) in prerequisites_data.dc_Transmission_nodes])
    else
        error("Unidentified grid model setting ($simulation_settings.dc_grid_model)")
        return -1
    end

    # converter variables
    JuMP.@variable(model, p_conv[g in collect(union(Set(prerequisites_data.ac_virtual_gen_set),Set(prerequisites_data.dc_virtual_gen_set)))])

    # Switching variables
    if simulation_settings.transmission_switching
        JuMP.@variable(model, z[l in prerequisites_data.switched_transmission_set], Bin)
    end

    if simulation_settings.substation_switching
        JuMP.@variable(model, z_l[l in prerequisites_data.reconf_set], Bin)
        JuMP.@variable(model, z_c[c in prerequisites_data.Coupler_set], Bin)
    end
end

function single_period_nodal_balance_ac_node!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Nodes_set = prereq.Nodes_set
    Branch_nodes = prereq.Branch_nodes
    load_set = prereq.load_set
    Gen_set = prereq.Gen_set
    ac_virtual_gen_set = prereq.ac_virtual_gen_set

    if formulation == :ACOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set],
            sum(model[:pij][i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_conv][g] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
        
        # REACTIVE NODAL BALANCE
        JuMP.@constraint(model, Qnodal[i in Nodes_set],
            sum(model[:qij][i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:q][g] for g in Gen_set if grid.Generators[g].GenBus_ID == i)
                 - sum(grid.Loads[l].Qd for l in load_set if grid.Loads[l].LoadBus_ID == i))
    elseif formulation == :DCOPF
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal[i in Nodes_set],
            sum(model[:pij][i,j] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
                - sum(grid.Loads[l].Pd for l in load_set if grid.Loads[l].LoadBus_ID == i) + sum(model[:p_conv][g] for g in ac_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_powerflow_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol; ideal_branches=[])
    
    Nodes_set = prereq.Nodes_set
    Branch_nodes = prereq.Branch_nodes
    unswitched_Transmission_nodes = prereq.unswitched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    Sbase = prereq.Sbase
    B = prereq.B

    if formulation == :ACOPF
        # ACTIVE POWER THROUGH LINE i-j
        JuMP.@NLconstraint(model, p_line[i in Nodes_set, j in Nodes_set; Set([i,j]) in unswitched_Transmission_nodes && branch_dictionary[Set([i,j])] ∉ ideal_branches],
            model[:pij][i,j] ==  Sbase*(model[:v][i]*model[:v][j]*(G[i,j]*cos(model[:δ][i]-model[:δ][j])+B[i,j]*sin(model[:δ][i]-model[:δ][j])) -(model[:v][i]^2)*G[i,j]))

        # REACTIVE POWER THROUGH LINE i-j
        JuMP.@NLconstraint(model, q_line[i in Nodes_set, j in Nodes_set; Set([i,j]) in unswitched_Transmission_nodes && branch_dictionary[Set([i,j])] ∉ ideal_branches],
            model[:pij][i,j] ==  Sbase*(model[:v][i]*model[:v][j]*(G[i,j]*sin(model[:δ][i]-model[:δ][j]) - B[i,j]*cos(model[:δ][i]-model[:δ][j])) +(model[:v][i]^2)*(B[i,j] - grid.Branches[branch_dictionary[Set([i,j])]].b/2)))

    elseif formulation == :DCOPF
        # ACTIVE POWER THROUGH LINE i-j
        JuMP.@constraint(model,pl[i in Nodes_set,j in Nodes_set; Set([i,j]) in unswitched_Transmission_nodes && branch_dictionary[Set([i,j])] ∉ ideal_branches],
            model[:pij][i,j] == Sbase*(B[i,j])*(model[:δ][i]-model[:δ][j]))
        
        # OPPOSITE FLOW CONSISTENCY
        JuMP.@constraint(model, pl_consistency[i in Nodes_set,j in Nodes_set; Set([i,j]) in Branch_nodes], 
        model[:pij][i,j] == -model[:pij][j,i])
    else
        error("Requested formulation ($formulation) is not implemented.")
    end  
end

function single_period_transmission_capacity_limits_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)

    Nodes_set = prereq.Nodes_set
    Sbase = prereq.Sbase
    branch_dictionary = prereq.branch_dictionary
    unswitched_Transmission_nodes = prereq.unswitched_Transmission_nodes

    if formulation == :ACOPF
        # LINE CAPACITY
        JuMP.@NLconstraint(model, Smax[i in Nodes_set, j in Nodes_set; Set([i,j]) in unswitched_Transmission_nodes],
            (model[:pij][i,j])^2 + (model[:qij][i,j])^2 ≤ ((Sbase)^2)*(grid.Branches[branch_dictionary[Set([i,j])]].rating)^2)
    elseif formulation == :DCOPF
        # LINE CAPACITY
        JuMP.@constraint(model,pl_rate[i in Nodes_set,j in Nodes_set; Set([i,j]) in unswitched_Transmission_nodes],
            -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating ≤ model[:pij][i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating)
    else
        error("Requested formulation ($formulation) is not implemented.")
    end  
    
end

function single_period_switched_transmission_capacity_limits_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model)
    
    Sbase = prereq.Sbase
    Nodes_set = prereq.Nodes_set
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    # LINE CAPACITY
    JuMP.@constraint(model,pl_rate_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*model[:z][branch_dictionary[Set([i,j])]] ≤ model[:pij][i,j] )

    JuMP.@constraint(model,pl_rate_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*model[:z][branch_dictionary[Set([i,j])]] )

end

function single_period_SG_switched_transmission_capacity_limits_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model,z)
    
    Sbase = prereq.Sbase
    Nodes_set = prereq.Nodes_set
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    branch_dictionary = prereq.branch_dictionary
    # LINE CAPACITY
    JuMP.@constraint(model,pl_rate_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])]] ≤ model[:pij][i,j] )

    JuMP.@constraint(model,pl_rate_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*z[branch_dictionary[Set([i,j])]] )

end

function single_period_switched_powerflow_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model; max_op = Inf)
    
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    Nodes_set = prereq.Nodes_set
    branch_dictionary = prereq.branch_dictionary
    B = prereq.B
    Sbase = prereq.Sbase
    M = grid.S_base*100
    JuMP.@constraint(model,pl_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] - Sbase*(B[i,j])*(model[:δ][i]-model[:δ][j]) ≤  (1-model[:z][branch_dictionary[Set([i,j])]])*M)
    
    JuMP.@constraint(model,pl_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] - Sbase*(B[i,j])*(model[:δ][i]-model[:δ][j]) ≥  -(1-model[:z][branch_dictionary[Set([i,j])]])*M)

    if max_op !== Inf
        JuMP.@constraint(model,maximum_allowed_switching_actions, sum(1-model[:z][branch_dictionary[l]] for l in switched_Transmission_nodes) ≤ max_op)
    end
end

function single_period_SG_switched_powerflow_ac_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, z)
    
    switched_Transmission_nodes = prereq.switched_Transmission_nodes
    Nodes_set = prereq.Nodes_set
    branch_dictionary = prereq.branch_dictionary
    B = prereq.B
    Sbase = prereq.Sbase
    M = grid.S_base*100
    JuMP.@constraint(model,pl_1[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] - Sbase*(B[i,j])*(model[:δ][i]-model[:δ][j]) ≤  (1-z[branch_dictionary[Set([i,j])]])*M)
    
    JuMP.@constraint(model,pl_2[i in Nodes_set,j in Nodes_set; Set([i,j]) in switched_Transmission_nodes],
        model[:pij][i,j] - Sbase*(B[i,j])*(model[:δ][i]-model[:δ][j]) ≥  -(1-z[branch_dictionary[Set([i,j])]])*M)
end

function single_period_reconf_split_constraints_ac_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model; max_reconf=Inf, max_splitting=Inf)
    
    aux_bus_set = prereq.aux_bus_set
    Reconf_dict = prereq.Reconf_dict
    Nodes_set = prereq.Nodes_set
    Reconf_nodes = prereq.Reconf_nodes
    Coupler_nodes = prereq.Coupler_nodes
    Coupler_dict = prereq.Coupler_dict
    Sbase = prereq.Sbase
    default_off_reconf = prereq.default_off_reconf
    Coupler_set = prereq.Coupler_set

    # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
    JuMP.@constraint(model, no_circular_path_constraint[bus in aux_bus_set],sum(model[:z_l][l] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
    M_δ = 2*π

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        model[:δ][j]-M_δ*(1-model[:z_l][Reconf_dict[ Set([i,j])]]) ≤ model[:δ][i] ) 

    JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes],
        model[:δ][i] ≤ model[:δ][j] + M_δ*(1-model[:z_l][Reconf_dict[ Set([i,j])]])) 

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        model[:δ][j] - M_δ*(1-model[:z_c][Coupler_dict[ Set([i,j])]]) ≤ model[:δ][i] ) 

    JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes],
        model[:δ][i] ≤ model[:δ][j] + M_δ*(1-model[:z_c][Coupler_dict[ Set([i,j])]]) ) 

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(model,reconf_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        -model[:z_l][Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating ≤ model[:pij][i,j])

    JuMP.@constraint(model,reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Reconf_nodes], 
        model[:pij][i,j] ≤ model[:z_l][Reconf_dict[Set([i,j])]] * Sbase * grid.Branches[Reconf_dict[Set([i,j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(model,coupler_cap_1[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        -model[:z_c][Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating ≤ model[:pij][i,j])

    JuMP.@constraint(model,coupler_cap_2[i in Nodes_set, j in Nodes_set; Set([i,j]) in Coupler_nodes], 
        model[:pij][i,j] ≤ model[:z_c][Coupler_dict[Set([i,j])]] * Sbase * grid.Branches[Coupler_dict[Set([i,j])]].rating)

    if max_reconf !== Inf
        JuMP.@constraint(model, max_reconf_actions, sum(model[:z_l][r] for r in default_off_reconf) ≤ max_reconf)
    end

    if max_splitting !== Inf
        JuMP.@constraint(model, max_splitting_actions, sum(model[:z_c][c] for c in Coupler_set) ≤ max_splitting)
    end

end

function single_period_reconf_split_SG_constraints_ac_grid!(prereq::OPF_Prerequisites, grid::PowerGrid, model::Model,z_l,z_c)

    aux_bus_set = prereq.aux_bus_set
    Reconf_dict = prereq.Reconf_dict
    Nodes_set = prereq.Nodes_set
    Reconf_nodes = prereq.Reconf_nodes
    Coupler_nodes = prereq.Coupler_nodes
    Coupler_dict = prereq.Coupler_dict
    Sbase = prereq.Sbase
    default_off_reconf = prereq.default_off_reconf
    Coupler_set = prereq.Coupler_set

    M_δ = 2 * π

    # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        model[:δ][j] - M_δ * (1 - z_l[Reconf_dict[Set([i, j])]]) ≤ model[:δ][i])

    JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        model[:δ][i] ≤ model[:δ][j] + M_δ * (1 - z_l[Reconf_dict[Set([i, j])]]))

    # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
    JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set; Set([i, j]) in Coupler_nodes],
        model[:δ][j] - M_δ * (1 - z_c[Coupler_dict[Set([i, j])]]) ≤ model[:δ][i])

    JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Coupler_nodes],
        model[:δ][i] ≤ model[:δ][j] + M_δ * (1 - z_c[Coupler_dict[Set([i, j])]]))

    # 2.7.1 Reconfiguration line capacity
    JuMP.@constraint(model, reconf_cap_1[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        -z_l[Reconf_dict[Set([i, j])]] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating ≤ model[:pij][i, j])

    JuMP.@constraint(model, reconf_cap_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Reconf_nodes],
        model[:pij][i, j] ≤ z_l[Reconf_dict[Set([i, j])]] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating)

    # 2.7.2 Reconfiguration line capacity
    JuMP.@constraint(model, coupler_cap_1[i in Nodes_set, j in Nodes_set; Set([i, j]) in Coupler_nodes],
        -z_c[Coupler_dict[Set([i, j])]] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating ≤ model[:pij][i, j])

    JuMP.@constraint(model, coupler_cap_2[i in Nodes_set, j in Nodes_set; Set([i, j]) in Coupler_nodes],
        model[:pij][i, j] ≤ z_c[Coupler_dict[Set([i, j])]] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating)

end

function single_period_generator_limits_ac_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Gen_set = prereq.Gen_set

    if formulation == :ACOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set], grid.Generators[g].Pg_min ≤ model[:p][g] ≤ grid.Generators[g].Pg_max)
        JuMP.@constraint(model, gen_reactive_power_limits[g in Gen_set], grid.Generators[g].Qg_min ≤ model[:q][g] ≤ grid.Generators[g].Qg_max)
    elseif formulation == :DCOPF
        # GENERATOR CAPACITY
        JuMP.@constraint(model, gen_active_power_limits[g in Gen_set], grid.Generators[g].Pg_min ≤ model[:p][g] ≤ grid.Generators[g].Pg_max)
        JuMP.@constraint(model, total_balance ,sum(model[:p][g] for g in Gen_set) == sum(grid.Loads[l].Pd for l in prereq.load_set))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_voltage_limits_ac_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    if formulation == :ACOPF
        Nodes_set = prereq.Nodes_set
        JuMP.@constraint(model, voltage_limits[i in Nodes_set], grid.Buses[i].V_min ≤ model[:v][i] ≤ grid.Buses[i].V_max)
    end
end

function single_period_angle_limits_ac_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model)
    Nodes_set = prereq.Nodes_set
    JuMP.@constraint(model, angle_limits[i in Nodes_set], grid.Buses[i].δ_min ≤ model[:δ][i] ≤ grid.Buses[i].δ_max)
    if !isnothing(prereq.reference_node)
        JuMP.@constraint(model,reference_node, model[:δ][prereq.reference_node] == 0)
    end
end

function single_period_converter_constraints!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    converter_set = prereq.converter_set
    b2b_gen_set = prereq.b2b_gen_set
    b2b_coupler_dict = prereq.b2b_coupler_dict
    ac_virtual_gen_set = prereq.ac_virtual_gen_set
    dc_virtual_gen_set = prereq.dc_virtual_gen_set
    
    if formulation == :Nonlinear
        error("NonLinear formulation for converters is not implemented yet")
    elseif formulation == :Linear_lossless
        JuMP.@constraint(model, conv_power[c in converter_set], model[:p_conv][grid.Converters[c].gen_ac_id] + model[:p_conv][grid.Converters[c].gen_dc_id] == 0)
        JuMP.@constraint(model,dc_link_power[link_id in keys(grid.DCLinks)], model[:p_conv][grid.DCLinks[link_id].Fr_gen_ID] + model[:p_conv][grid.DCLinks[link_id].To_gen_ID] == 0)
        JuMP.@constraint(model, conv_cap_ac_side[g in ac_virtual_gen_set], grid.Generators[g].Pg_min ≤ model[:p_conv][g] ≤ grid.Generators[g].Pg_max)
        JuMP.@constraint(model, conv_cap_dc_side[g in dc_virtual_gen_set], grid.Generators[g].Pg_min ≤ model[:p_conv][g] ≤ grid.Generators[g].Pg_max)
        if b2b_gen_set !== []
            # HYBRID SPLITTING CONSTRAINTS
            JuMP.@constraint(model,VSC_cap1[g in b2b_gen_set],(1-model[:z_c][b2b_coupler_dict[g]])*grid.Generators[g].Pg_min ≤ model[:p_conv][g])
            JuMP.@constraint(model,VSC_cap2[g in b2b_gen_set], model[:p_conv][g] ≤ (1-model[:z_c][b2b_coupler_dict[g]])*grid.Generators[g].Pg_max)
        end

    elseif formulation == :Linear_lossy
        error("Linear_lossy formulation for converters is not implemented yet")
        # JuMP.@constraint(model, conv_power[c in converter_set], p_conv[grid.Converters[c].gen_ac_id] + p_conv[grid.Converters[c].gen_dc_id] == losses)
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_nodal_balance_dc_node!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes
    dc_virtual_gen_set = prereq.dc_virtual_gen_set
    
    if formulation == :NonLinear
        error("NonLinear nodal balance is not implemented yet for DC grids")
    elseif formulation == :Linear
        # ACTIVE NODAL BALANCE
        JuMP.@constraint(model, Pnodal_dc[i in dc_Nodes_set],
            sum(model[:pij_dc][i,j] for j = dc_Nodes_set if Set([i,j]) in dc_Transmission_nodes) == sum(model[:p_conv][g] for g in dc_virtual_gen_set if grid.Generators[g].GenBus_ID == i))
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function single_period_transmission_capacity_limits_dc_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    Sbase = prereq.Sbase
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes
    dc_branch_dictionary = prereq.dc_branch_dictionary

    if formulation == :Linear
        # LINE CAPACITY
        JuMP.@constraint(model,pl_rate_dc[i in dc_Nodes_set,j in dc_Nodes_set; Set([i,j]) in dc_Transmission_nodes],
            -Sbase*grid.DCBranches[dc_branch_dictionary[Set([i,j])]].rating ≤ model[:pij_dc][i,j] ≤ Sbase*grid.DCBranches[dc_branch_dictionary[Set([i,j])]].rating)
    elseif formulation == :NonLinear
        error("NonLinear flow limits in DC grid is not implemented yet")
    else
        error("Requested formulation ($formulation) is not implemented")
    end
end

function single_period_powerflow_dc_branch!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    
    dc_Nodes_set = prereq.dc_Nodes_set
    dc_Transmission_nodes = prereq.dc_Transmission_nodes

    if formulation == :Linear
        JuMP.@constraint(model,dc_Pij[i in dc_Nodes_set, j in dc_Nodes_set;  Set([i,j]) in dc_Transmission_nodes], model[:pij_dc][i,j] + model[:pij_dc][j,i] == 0 )
    elseif formulation == :NonLinear
        error("NonLinear power flow formulation for DC grid is not implemented yet")
    else
        error("Requested formulation ($formulation) is not implemented")
    end
end

function single_period_voltage_limits_dc_grid!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model, formulation ::Symbol)
    if formulation == :NonLinear
        dc_Nodes_set = prereq.dc_Nodes_set
        JuMP.@constraint(model, dc_voltage[i in dc_Nodes_set], grid.DCBuses[i].V_min ≤ model[:v_dc][i] ≤ grid.DCBuses[i].V_max)
    end
end

function single_period_objective!(prereq ::OPF_Prerequisites, grid ::PowerGrid, model ::Model,transmission_switching ::Bool, substation_switching ::Bool)
    
    if transmission_switching && substation_switching
        Coupler_set = prereq.Coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g]+grid.Generators[g].C0 for g in prereq.Gen_set)+sum(1-model[:z_c][c] for c in Coupler_set)*1e-4+sum(1-model[:z][l] for l in prereq.switched_transmission_set)*1e-4)
    elseif transmission_switching && ! substation_switching
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g]+grid.Generators[g].C0 for g in prereq.Gen_set)+sum(1-model[:z][l] for l in prereq.switched_transmission_set)*1e-4)
    elseif substation_switching && ! transmission_switching
        Coupler_set = prereq.Coupler_set
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g]+grid.Generators[g].C0 for g in prereq.Gen_set)+sum(1-model[:z_c][c] for c in Coupler_set)*1e-4)
    else
        JuMP.@objective(model,Min,sum(grid.Generators[g].C1*model[:p][g]+grid.Generators[g].C0 for g in prereq.Gen_set))
    end

end