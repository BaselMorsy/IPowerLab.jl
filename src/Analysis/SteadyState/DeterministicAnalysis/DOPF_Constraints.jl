@with_kw mutable struct DOPF_SimulationSettings
    time_horizon = [1] # can be any range e.g. 1:24, note that any non-continous range can pose problems for temporal transition constraints
    ac_grid_model = :Bθ # :AC
    dc_grid_model = :NF # NF stands for network flow and this is currently the only available model but soon :FNL (full nonlinear) model will be implemented
    converter_model = :NF_lossless # means that P_conv_ac + P_conv_dc = 0, can be any thing from [:NF_lossy, :FNL] but only :NF_lossles implemented now
    dynamic_converter_control = true # allows for converter redispatch post-contingency
    converter_modularization = :discrete # useful in case of busbar splitting, where converter modules can be modeled using binary or continous variables
    
    load_shedding = [:post] # [:pre], [:post], [:pre,:post] -> if empty, default behaviour is [:pre,:post]
    transmission_switching = [:post] # [:pre] , [:post] , [:pre,:post] -> if empty, default behaviour is [:pre,:post]
    substation_switching = Dict("splitting" => [:post], "reconf" => [:pre]) # -> if empty, default behaviour is [:pre,:post]
    activate_temporal_switching = false
    max_transmission_switching = Dict("pre_contingency" => Inf, "post_contingency" => Inf, "MCDT" => 1)
    max_substation_reconf = Dict("MCDT" => 1)
    max_busbar_splitting = Dict("pre_contingency" => Inf, "post_contingency" => Inf, "MCDT" => 1)
    
    contingency_types = [:ac_branch, :dc_branch, :dc_link, :ac_gen, :dc_gen, :conv, :coupler] # Can be any subset of the provided default vector
    contingency_redispatch_condition = :loss_of_generation # could also be :none or :all

    NLP_solver = [] #Ipopt.Optimizer
    MILP_solver = [] #Gurobi.Optimizer
    Meta_solver = :none # => A pipeline of solution algorithms. To be determined. Here is where innovation happens!
end

@with_kw mutable struct DOPF_Prerequisites
    time_horizon
    base_MVA
    ac_node_ids # all ac nodes ids
    ac_aux_bus_ids # auxilliary node ids
    ac_branch_ids # all branch ids including auxilliary branches
    ac_gen_ids # all real generator ids (virtual generators not included)
    ac_load_ids # load ids
    ac_load_shedding_ids # default -> same as ac_load_ids

    ac_gen_id_to_gen_root # dictionary g -> g_root where g_root is the main generator NOT the duplicate one in case of different down-reg costs
    root_gen_to_duplicate_gen # dictionary g_root -> g_duplicate
    contingency_redispatch # contingency indices i where (i > 1) at which generation is allowed to redispatch
    commitable_gen_ids # generators that can be commited (u_gt)
    non_commitable_gen_ids # generators that are deemed to run non-stop
    fixed_commitments # dictionary of generators that are already committed due to contracts or anything (u_gt = 1 or 0)
    fixed_schedules # dictionary of fixed schedules (p_gtk = x)

    ac_static_branch_ids # ids of unswitchable real branches (no auxilliary branches)
    ac_active_dynamic_branch_ids # ids of switchable real branches
    ac_active_reconf_ids # ids of active switchable auxilliary (reconfiguration) branches
    ac_active_coupler_ids # # ids of active busbar coupler auxilliary branches
    ac_fixed_dynamic_branch_ids # ids of real switchable but fixed status branches
    ac_fixed_reconf_ids # ids of auxilliary switchable but fixed branches
    ac_fixed_coupler_ids # ids of auxilliary siwtchable but fixed couplers
    fixed_topology # dictionary of fixed topology [line_id] -> status (1,0)
    normally_opened_reconf_lines # list of normally open reconf lines 

    dc_link_ids
    dc_link_vritual_gen_ids
    conv_ids
    b2b_conv_ids
    b2b_gen_ids
    conv_ac_side_virtual_gen_ids
    conv_dc_side_virtual_gen_ids
    ac_gen_to_conv_id_dict
    dc_gen_to_conv_id_dict
    conv_duplets # dictionary of converter duplets

    dc_node_ids
    dc_branch_ids
    dc_gen_ids
    dc_load_ids

    Schedule # comes from market, if non zero means that DOPF is a redispatch clearing
    Order_Book # comes from market, can be filled from generators' data in `grid.Generators` or can be filled arbitrarily outside (e.g. using a bidding agent)
    
    Contingency_Map
    k # vector of contingency indices
    k_t # dictionary of contingencies considered at time t

    M_l # Dictionary of M_l values, where the keys are the AC branch ids (Big-M for switchable branches)
    M_δ # scalar number for the maximum Δδ between any two sections in a substation (Big-M for reconfiguration and splitting)
    M_E # scalar number for auxilliary branches capacities 

    relaxed_physics_lines # lines that doesn't comply to power flow equations and can have arbitrary flow as long as it doesn't violate their capacity limits
    relaxed_physics_nodes # nodes that has a zero phase angle and 1 voltage magnitude   
    relaxed_capacity_lines # lines that has no capacity constraints

    reference_node = nothing
end

function DOPF_variable_initialization!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    """
    Variable initialization for DOPF simulation, based on simulation settings `simulation_settings` and prerequisites `prerequisites_data`
    """
    sides = [1,2] # two sides of the power flow (fr,to) and (to,fr)
    # variable initialization
    if simulation_settings.ac_grid_model == :Bθ

        ###################################################################################
        JuMP.@variable(model, δ[i in prerequisites_data.ac_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_branch_ac[l in prerequisites_data.ac_branch_ids, s in sides, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])

        ###################################################################################
        # Binary commitment variables
        if length(prerequisites_data.commitable_gen_ids) != 0
            JuMP.@variable(model, u_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin)
            JuMP.@variable(model, α_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # start-up
            JuMP.@variable(model, β_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # shut-down
        end

        # Fixed commitment variables
        if length(keys(prerequisites_data.fixed_commitments)) != 0
            JuMP.@variable(model, u_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        # Branch switching binary variables
        if prerequisites_data.ac_active_dynamic_branch_ids != []
            JuMP.@variable(model, z_l[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]], Bin)
            if simulation_settings.activate_temporal_switching
                JuMP.@variable(model, α_l[l in prerequisites_data.ac_active_dynamic_branch_ids, t in prerequisites_data.time_horizon], Bin)
                JuMP.@variable(model, β_l[l in prerequisites_data.ac_active_dynamic_branch_ids, t in prerequisites_data.time_horizon], Bin)
            end
        end
        # Branch switching fixed continous variables
        if prerequisites_data.ac_fixed_dynamic_branch_ids != []
            JuMP.@variable(model, z_l_f[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        end
        ###################################################################################
        # Safe guards for substation switching
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
        ###################################################################################
        # Substation switching binary variables
        ###################################################################################
        if prerequisites_data.ac_active_reconf_ids != []
            JuMP.@variable(model, z_r[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]], Bin)
            if simulation_settings.activate_temporal_switching
                JuMP.@variable(model, α_r[l in prerequisites_data.ac_active_reconf_ids, t in prerequisites_data.time_horizon], Bin)
                JuMP.@variable(model, β_r[l in prerequisites_data.ac_active_reconf_ids, t in prerequisites_data.time_horizon], Bin)
            end
        end

        if prerequisites_data.ac_active_coupler_ids != []
            JuMP.@variable(model, z_c[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]], Bin)
            if simulation_settings.activate_temporal_switching
                JuMP.@variable(model, α_c[l in prerequisites_data.ac_active_coupler_ids, t in prerequisites_data.time_horizon], Bin)
                JuMP.@variable(model, β_c[l in prerequisites_data.ac_active_coupler_ids, t in prerequisites_data.time_horizon], Bin)
            end
        end
        ###################################################################################
        if prerequisites_data.ac_fixed_reconf_ids != []
            JuMP.@variable(model, z_r_f[r in prerequisites_data.ac_fixed_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        end

        if prerequisites_data.ac_fixed_coupler_ids != []
            JuMP.@variable(model, z_c_f[c in prerequisites_data.ac_fixed_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        end
        ###################################################################################

    elseif simulation_settings.ac_grid_model == :AC

        JuMP.@variable(model, v[i in prerequisites_data.ac_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, δ[i in prerequisites_data.ac_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, p_branch_ac[l in prerequisites_data.ac_branch_ids, s in sides, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, q_branch_ac[l in prerequisites_data.ac_branch_ids, s in sides, k in prerequisites_data.k, t in prerequisites_data.time_horizon])

        JuMP.@variable(model, p_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, q_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        
        if length(prerequisites_data.ac_load_shedding_ids) != 0
            JuMP.@variable(model, p_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
            JuMP.@variable(model, q_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        # Binary commitment variables
        if length(prerequisites_data.commitable_gen_ids) != 0
            JuMP.@variable(model, u_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon])
            JuMP.@variable(model, α_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon]) # start-up
            JuMP.@variable(model, β_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon]) # shut-down
        end

        # Fixed commitment variables
        if length(keys(prerequisites_data.fixed_commitments)) != 0
            JuMP.@variable(model, u_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        # Branch switching binary variables
        if prerequisites_data.ac_active_dynamic_branch_ids != []
            JuMP.@variable(model, z_l[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        # Branch switching fixed continous variables
        if prerequisites_data.ac_fixed_dynamic_branch_ids != []
            JuMP.@variable(model, z_l_f[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        # Safe guards for substation switching
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
        ###################################################################################
        # Substation switching binary variables
        ###################################################################################
        if prerequisites_data.ac_active_reconf_ids != []
            JuMP.@variable(model, z_r[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end

        if prerequisites_data.ac_active_coupler_ids != []
            JuMP.@variable(model, z_c[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        if prerequisites_data.ac_fixed_reconf_ids != []
            JuMP.@variable(model, z_r_f[r in prerequisites_data.ac_fixed_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end

        if prerequisites_data.ac_fixed_coupler_ids != []
            JuMP.@variable(model, z_c_f[c in prerequisites_data.ac_fixed_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        ###################################################################################
    else
        error("Invalid AC grid model setting ($(simulation_settings.ac_grid_model))")
        return -1
    end

    # Sanity check
    if length(prerequisites_data.dc_branch_ids) != 0 && length(prerequisites_data.dc_node_ids) == 0
        error("Cannot have DC branches with no DC nodes!")
    elseif length(prerequisites_data.dc_gen_ids) != 0 && length(prerequisites_data.dc_node_ids) == 0
        error("Cannot have DC generators with no DC nodes!")
    elseif length(prerequisites_data.dc_load_ids) != 0 && length(prerequisites_data.dc_node_ids) == 0
        error("Cannot have DC loads with no DC nodes!")
    end

    if simulation_settings.dc_grid_model == :NF
        JuMP.@variable(model, p_branch_dc[l in prerequisites_data.dc_branch_ids, s in sides, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_gen_dc[g in prerequisites_data.dc_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
    elseif simulation_settings.dc_grid_model == :FNL
        JuMP.@variable(model, v_dc[i in prerequisites_data.dc_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, p_branch_dc[l in prerequisites_data.dc_branch_ids, s in sides, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, p_gen_dc[g in prerequisites_data.dc_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
    else
        error("Invalid DC grid model setting ($(simulation_settings.dc_grid_model))")
        return -1
    end
    # maybe we don't need all the if statements down there !!
    if length(prerequisites_data.conv_dc_side_virtual_gen_ids) != 0
        JuMP.@variable(model, p_conv_dc[conv_id in prerequisites_data.conv_dc_side_virtual_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_conv_ac[conv_id in prerequisites_data.conv_ac_side_virtual_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
    end
    if length(prerequisites_data.b2b_gen_ids) != 0
        JuMP.@variable(model, p_conv_b2b[conv_id in prerequisites_data.b2b_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
    end
    if simulation_settings.converter_modularization == :continuous && length(collect(keys(prerequisites_data.conv_duplets))) != 0
        JuMP.@variable(model, 0 ≤ γ_conv[duplet_id in collect(keys(prerequisites_data.conv_duplets)), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]] ≤ 1)
    end

    if length(prerequisites_data.dc_link_ids) != 0
        # DC link variables
        JuMP.@variable(model, p_dclink[link_id in prerequisites_data.dc_link_vritual_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
    end

end

function DOPF_nodal_balance_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        
        JuMP.@constraint(model, ac_active_nodal_balance[n in prerequisites_data.ac_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                sum(model[:p_gen_ac][g,k,t] + get(get(prerequisites_data.Schedule["gen"],g,Dict()),t,0) for g in grid.Buses[n].ConnectedGensIDs if g in prerequisites_data.ac_gen_ids)
                + sum([model[:p_conv_b2b][g,k,t] for g in grid.Buses[n].ConnectedGensIDs if g in prerequisites_data.b2b_gen_ids], init = 0)
                + sum([model[:p_conv_ac][g,k,t] for g in grid.Buses[n].ConnectedGensIDs if g in prerequisites_data.conv_ac_side_virtual_gen_ids], init = 0)
                + sum([model[:p_dc_link][g,k,t] for g in grid.Buses[n].ConnectedGensIDs if g in prerequisites_data.dc_link_vritual_gen_ids], init = 0)
                + sum([model[:p_ls_ac][d,k,t] for d in grid.Buses[n].ConnectedLoadsIDs if d in prerequisites_data.ac_load_shedding_ids], init = 0)
                - sum([grid.Loads[d].Pd_t[t] for d in grid.Buses[n].ConnectedLoadsIDs if d in prerequisites_data.ac_load_ids], init = 0)
                == sum([model[:p_branch_ac][l,1,k,t] for l in grid.Buses[n].ConnectedLinesIDs if grid.Branches[l].Fr_bus_ID == n], init = 0)
                + sum([model[:p_branch_ac][l,2,k,t] for l in grid.Buses[n].ConnectedLinesIDs if grid.Branches[l].To_bus_ID == n], init = 0)
                )
    elseif simulation_settings.ac_grid_model == :AC
        if length(prerequisites_data.ac_load_shedding_ids) != 0
        else
        end
    end
end

function DOPF_powerflow_ac_branch_static!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, ac_active_flow_fr_to[l in prerequisites_data.ac_static_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,1,k,t] == prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].Fr_bus_ID,k,t] - model[:δ][grid.Branches[l].To_bus_ID,k,t]))
        
        JuMP.@constraint(model, ac_active_flow_to_fr[l in prerequisites_data.ac_static_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,2,k,t] == prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].To_bus_ID,k,t] - model[:δ][grid.Branches[l].Fr_bus_ID,k,t]))
        
        if length(prerequisites_data.relaxed_physics_lines) != 0
            JuMP.@constraint(model, ac_branch_relaxed_consistent_flows[l in prerequisites_data.relaxed_physics_lines, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_branch_ac][l,1,k,t] + model[:p_branch_ac][l,2,k,t] == 0)
        end
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_thermal_limits_ac_branch_static!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, ac_active_flow_limits[l in prerequisites_data.ac_static_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_capacity_lines && k ∈ prerequisites_data.k_t[t]],
            -grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA ≤ model[:p_branch_ac][l,s,k,t] ≤ grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA)
    elseif simulation_settings.ac_grid_model == :AC
        error("AC model is not implemented yet")
    else
        error("AC grid model $(simulation_settings.ac_grid_mode) is invalid.")
    end
end

function DOPF_powerflow_ac_branch_dynamic!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        # upper limit (fr-to)
        JuMP.@constraint(model, ac_switched_active_flow_fr_to_up[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].Fr_bus_ID,k,t] - model[:δ][grid.Branches[l].To_bus_ID,k,t]) - model[:p_branch_ac][l,1,k,t] ≤ (1-model[:z_l][l,k,t])*prerequisites_data.M_l[l] )
        # lower limit (fr-to)
        JuMP.@constraint(model, ac_switched_active_flow_fr_to_down[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].Fr_bus_ID,k,t] - model[:δ][grid.Branches[l].To_bus_ID,k,t]) - model[:p_branch_ac][l,1,k,t] ≥ -(1-model[:z_l][l,k,t])*prerequisites_data.M_l[l])
        
        # upper limit (to-fr)
        JuMP.@constraint(model, ac_switched_active_flow_to_fr_up[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].To_bus_ID,k,t] - model[:δ][grid.Branches[l].Fr_bus_ID,k,t]) - model[:p_branch_ac][l,2,k,t] ≤ (1-model[:z_l][l,k,t])*prerequisites_data.M_l[l])
        
        # lower limit (to-fr)
        JuMP.@constraint(model, ac_switched_active_flow_to_fr_down[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].To_bus_ID,k,t] - model[:δ][grid.Branches[l].Fr_bus_ID,k,t]) - model[:p_branch_ac][l,2,k,t] ≥ -(1-model[:z_l][l,k,t])*prerequisites_data.M_l[l])
            
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_thermal_limits_ac_branch_dynamic!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, ac_switchable_active_flow_limits_up[l in prerequisites_data.ac_active_dynamic_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,s,k,t] ≤ model[:z_l][l,k,t]*grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA)
        
        JuMP.@constraint(model, ac_switchable_active_flow_limits_down[l in prerequisites_data.ac_active_dynamic_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,s,k,t] ≥ -model[:z_l][l,k,t]*grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA)
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_powerflow_ac_branch_dynamic_fixed!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        # upper limit (fr-to)
        JuMP.@constraint(model, ac_switched_active_flow_fr_to_up_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].Fr_bus_ID,k,t] - model[:δ][grid.Branches[l].To_bus_ID,k,t]) - model[:p_branch_ac][l,1,k,t] ≤ (1-model[:z_l_f][l,k,t])*prerequisites_data.M_l[l] )
        # lower limit (fr-to)
        JuMP.@constraint(model, ac_switched_active_flow_fr_to_down_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].Fr_bus_ID,k,t] - model[:δ][grid.Branches[l].To_bus_ID,k,t]) - model[:p_branch_ac][l,1,k,t] ≥ -(1-model[:z_l_f][l,k,t])*prerequisites_data.M_l[l])
        
        # upper limit (to-fr)
        JuMP.@constraint(model, ac_switched_active_flow_to_fr_up_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].To_bus_ID,k,t] - model[:δ][grid.Branches[l].Fr_bus_ID,k,t]) - model[:p_branch_ac][l,2,k,t] ≤ (1-model[:z_l_f][l,k,t])*prerequisites_data.M_l[l])
        
        # lower limit (to-fr)
        JuMP.@constraint(model, ac_switched_active_flow_to_fr_down_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_physics_lines && k ∈ prerequisites_data.k_t[t]],
            prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA*(1/grid.Branches[l].x)*(model[:δ][grid.Branches[l].To_bus_ID,k,t] - model[:δ][grid.Branches[l].Fr_bus_ID,k,t]) - model[:p_branch_ac][l,2,k,t] ≥ -(1-model[:z_l_f][l,k,t])*prerequisites_data.M_l[l])
            
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_thermal_limits_ac_branch_dynamic_fixed!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, ac_switchable_active_flow_limits_up_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_capacity_lines && k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,s,k,t] ≤ model[:z_l_f][l,k,t]*grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA)
        
        JuMP.@constraint(model, ac_switchable_active_flow_limits_down_fixed[l in prerequisites_data.ac_fixed_dynamic_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; l ∉ prerequisites_data.relaxed_capacity_lines && k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][l,s,k,t] ≥ -model[:z_l_f][l,k,t]*grid.Branches[l].rating*prerequisites_data.Contingency_Map["ac_branch"][l,k]*prerequisites_data.base_MVA)
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_substation_switching_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ

        # Switching constraint to avoid connecting an element to two busbar sections at the same time
        JuMP.@constraint(model, no_circular_path_constraint[bus in prerequisites_data.ac_aux_bus_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            sum(model[:z_r][l, k, t] for l in grid.Buses[bus].ConnectedLinesIDs if (grid.Branches[l].BranchType == 1) && (l in prerequisites_data.ac_active_reconf_ids)) == 1)
        #######################################################################################################################
        # Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
        JuMP.@constraint(model, phase_angle_equivalence_reconf_up[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[r].Fr_bus_ID, k, t] - model[:δ][grid.Branches[r].To_bus_ID, k, t] ≤ prerequisites_data.M_δ * (1 - model[:z_r][r, k, t]))

        JuMP.@constraint(model, phase_angle_equivalence_reconf_down[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[r].Fr_bus_ID, k, t] - model[:δ][grid.Branches[r].To_bus_ID, k, t] ≥ -prerequisites_data.M_δ * (1 - model[:z_r][r, k, t]))
        #######################################################################################################################
        # Phase angle constraints accross all couplers to be the same if the switch is 1
        JuMP.@constraint(model, phase_angle_equivalence_coupler_up[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[c].Fr_bus_ID, k, t] - model[:δ][grid.Branches[c].To_bus_ID, k, t] ≤ prerequisites_data.M_δ * (1 - model[:z_c][c, k, t]*prerequisites_data.Contingency_Map["coupler"][c,k]))

        JuMP.@constraint(model, phase_angle_equivalence_coupler_down[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[c].Fr_bus_ID, k, t] - model[:δ][grid.Branches[c].To_bus_ID, k, t] ≥ -prerequisites_data.M_δ * (1 - model[:z_c][c, k, t]*prerequisites_data.Contingency_Map["coupler"][c,k]))
        #######################################################################################################################
        # Reconfiguration line capacity
        JuMP.@constraint(model, reconf_flow[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][r, 1, k, t] + model[:p_branch_ac][r, 2, k, t] == 0)

        JuMP.@constraint(model, reconf_cap_up[r in prerequisites_data.ac_active_reconf_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][r, s, k, t] ≤ model[:z_r][r, k, t] * prerequisites_data.M_E)
        JuMP.@constraint(model, reconf_cap_down[r in prerequisites_data.ac_active_reconf_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][r, s, k, t] ≥ -model[:z_r][r, k, t] * prerequisites_data.M_E)
        #######################################################################################################################
        # Couplers capacity
        JuMP.@constraint(model, coupler_flow[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][c, 1, k, t] + model[:p_branch_ac][c, 2, k, t] == 0)

        JuMP.@constraint(model, coupler_cap_up[c in prerequisites_data.ac_active_coupler_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][c, s, k, t] ≤ model[:z_c][c, k, t] * prerequisites_data.M_E * prerequisites_data.Contingency_Map["coupler"][c,k])
        JuMP.@constraint(model, coupler_cap_down[c in prerequisites_data.ac_active_coupler_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][c, s, k, t] ≥ -model[:z_c][c, k, t] * prerequisites_data.M_E * prerequisites_data.Contingency_Map["coupler"][c,k])
    elseif simulation_settings.ac_grid_model == :AC
    end

    # Symmetry reduction constraints
    if length(keys(grid.Substations)) != 0
        JuMP.@constraint(model, symmetry_reduction[s in collect(keys(grid.Substations)), k in prerequisites_data.k, t in prerequisites_data.time_horizon; length(grid.Substations[s].BusbarSections_IDs) == 2 && k ∈ prerequisites_data.k_t[t]],
            sum(model[:z_r][r,k,t] for r in grid.Buses[grid.Substations[s].BusbarSections_IDs[1]].ConnectedLinesIDs if grid.Branches[r].BranchType == 1) ≥ sum(model[:z_r][r,k,t] for r in grid.Buses[grid.Substations[s].BusbarSections_IDs[2]].ConnectedLinesIDs if grid.Branches[r].BranchType == 1))
    end
end

function DOPF_substation_switching_ac_node_fixed!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ

        # Switching constraint to avoid connecting an element to two busbar sections at the same time
        JuMP.@constraint(model, no_circular_path_constraint_fixed[bus in prerequisites_data.ac_aux_bus_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            sum(model[:z_r_f][l, k, t] for l in grid.Buses[bus].ConnectedLinesIDs if (grid.Branches[l].BranchType == 1) && (l in prerequisites_data.ac_active_reconf_ids)) == 1)
        #######################################################################################################################
        # Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
        JuMP.@constraint(model, phase_angle_equivalence_reconf_up_fixed[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[r].Fr_bus_ID, k, t] - model[:δ][grid.Branches[r].To_bus_ID, k, t] ≤ prerequisites_data.M_δ * (1 - model[:z_r_f][r, k, t]))

        JuMP.@constraint(model, phase_angle_equivalence_reconf_down_fixed[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[r].Fr_bus_ID, k, t] - model[:δ][grid.Branches[r].To_bus_ID, k, t] ≥ -prerequisites_data.M_δ * (1 - model[:z_r_f][r, k, t]))
        #######################################################################################################################
        # Phase angle constraints accross all couplers to be the same if the switch is 1
        JuMP.@constraint(model, phase_angle_equivalence_coupler_up_fixed[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[c].Fr_bus_ID, k, t] - model[:δ][grid.Branches[c].To_bus_ID, k, t] ≤ prerequisites_data.M_δ * (1 - model[:z_c_f][c, k, t]*prerequisites_data.Contingency_Map["coupler"][c,k]))

        JuMP.@constraint(model, phase_angle_equivalence_coupler_down_fixed[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:δ][grid.Branches[c].Fr_bus_ID, k, t] - model[:δ][grid.Branches[c].To_bus_ID, k, t] ≥ -prerequisites_data.M_δ * (1 - model[:z_c_f][c, k, t]*prerequisites_data.Contingency_Map["coupler"][c,k]))
        #######################################################################################################################
        # Reconfiguration line capacity
        JuMP.@constraint(model, reconf_cap_up_fixed[r in prerequisites_data.ac_active_reconf_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][r, s, k, t] ≤ model[:z_r_f][r, k, t] * prerequisites_data.M_E)
        JuMP.@constraint(model, reconf_cap_down_fixed[r in prerequisites_data.ac_active_reconf_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][r, s, k, t] ≥ -model[:z_r_f][r, k, t] * prerequisites_data.M_E)
        #######################################################################################################################
        # Couplers capacity
        JuMP.@constraint(model, coupler_cap_up_fixed[c in prerequisites_data.ac_active_coupler_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][c, s, k, t] ≤ model[:z_c_f][c, k, t] * prerequisites_data.M_E * prerequisites_data.Contingency_Map["coupler"][c,k])
        JuMP.@constraint(model, coupler_cap_down_fixed[r in prerequisites_data.ac_active_coupler_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_ac][c, s, k, t] ≥ -model[:z_c_f][c, k, t] * prerequisites_data.M_E * prerequisites_data.Contingency_Map["coupler"][c,k])
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_generator_limits_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        
        if length(prerequisites_data.non_commitable_gen_ids) != 0
            JuMP.@constraint(model, ac_gen_cap_non_commitable[g in prerequisites_data.non_commitable_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][2]*prerequisites_data.Contingency_Map["ac_gen"][g,k] ≤ model[:p_gen_ac][g,k,t] 
                ≤ prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][1]*prerequisites_data.Contingency_Map["ac_gen"][g,k])
        end

        if length(prerequisites_data.commitable_gen_ids) != 0
            JuMP.@constraint(model, ac_gen_cap_commitable_up[g in prerequisites_data.commitable_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,k,t] ≤ prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][1]*model[:u_gt][g,t]*prerequisites_data.Contingency_Map["ac_gen"][g,k])
            
            JuMP.@constraint(model, ac_gen_cap_commitable_down[g in prerequisites_data.commitable_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,k,t] ≥ prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][2]*model[:u_gt][g,t]*prerequisites_data.Contingency_Map["ac_gen"][g,k])
        end

        if length(keys(prerequisites_data.fixed_commitments)) != 0
            JuMP.@constraint(model, ac_gen_cap_fixed_commitments_up[g in keys(prerequisites_data.fixed_commitments), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,k,t] ≤ prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][1]*model[:u_gt_f][g,t]*prerequisites_data.Contingency_Map["ac_gen"][g,k])
                
            JuMP.@constraint(model, ac_gen_cap_fixed_commitments_down[g in keys(prerequisites_data.fixed_commitments), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,k,t] ≥ prerequisites_data.Order_Book.Gen_bids[g]["qty"][t][2]*model[:u_gt_f][g,t]*prerequisites_data.Contingency_Map["ac_gen"][g,k])
        end

    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_load_shedding_limits_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, load_shedding_limits_up[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t] && prerequisites_data.Order_Book.Load_bids[d]["qty"][t] ≥ 0 ],
            0 ≤ model[:p_ls_ac][d,k,t] ≤ prerequisites_data.Order_Book.Load_bids[d]["qty"][t])

        JuMP.@constraint(model, load_shedding_limits_down[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t] && prerequisites_data.Order_Book.Load_bids[d]["qty"][t] ≤ 0],
            prerequisites_data.Order_Book.Load_bids[d]["qty"][t] ≤ model[:p_ls_ac][d,k,t] ≤ 0)
    elseif simulation_settings.ac_grid_model == :AC
        if length(prerequisites_data.ac_load_shedding_ids) != 0
        else
        end
    end
end

function DOPF_voltage_limits_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.ac_grid_model == :Bθ
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_angle_limits_ac_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites) 
    if simulation_settings.ac_grid_model == :Bθ
        JuMP.@constraint(model, angle_limits[i in prerequisites_data.ac_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            grid.Buses[i].δ_min ≤ model[:δ][i,k,t] ≤ grid.Buses[i].δ_max)
        if length(prerequisites_data.relaxed_physics_nodes) != 0
            JuMP.@constraint(model, relaxed_nodes_angles[i in prerequisites_data.relaxed_physics_nodes, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:δ][i,k,t] == 0)
        end
        if !isnothing(prerequisites_data.reference_node)
            JuMP.@constraint(model, reference_node[k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:δ][prerequisites_data.reference_node,k,t] == 0)
        end
    elseif simulation_settings.ac_grid_model == :AC
    end
end

function DOPF_converter_flow!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.converter_model == :FNL
        error("NonLinear converter model is not implemented yet")
    elseif simulation_settings.converter_model == :NF_lossless

        JuMP.@constraint(model, converter_flow[c in prerequisites_data.conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_dc][grid.Converters[c].gen_dc_id, k, t] == 0)
        
        JuMP.@constraint(model, b2b_converter_flow[c in prerequisites_data.b2b_conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_conv_b2b][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_b2b][grid.Converters[c].gen_dc_id, k, t] == 0)
    elseif simulation_settings.converter_model == :NF_lossy
        
        JuMP.@constraint(model, converter_flow_1[c in prerequisites_data.conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_dc][grid.Converters[c].gen_dc_id, k, t] == grid.Converters[c].loss_a + grid.Converters[c].loss_b*model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t])

        # JuMP.@constraint(model, converter_flow_2[c in prerequisites_data.conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
        #     model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_dc][grid.Converters[c].gen_dc_id, k, t] ≥ grid.Converters[c].loss_a - grid.Converters[c].loss_b*model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t])

        # JuMP.@constraint(model, converter_flow_3[c in prerequisites_data.conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
        #     model[:p_conv_ac][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_dc][grid.Converters[c].gen_dc_id, k, t] ≤ grid.Converters[c].loss_a + grid.Converters[c].loss_b*grid.Converters[c].rate)
        
        JuMP.@constraint(model, b2b_converter_flow[c in prerequisites_data.b2b_conv_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_conv_b2b][grid.Converters[c].gen_ac_id, k, t] + model[:p_conv_b2b][grid.Converters[c].gen_dc_id, k, t] == 0)
    else
        error("Requested converter model ($(simulation_settings.converter_model)) is not implemented.")
    end
end

function DOPF_converter_capacity!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.converter_model == :FNL
        error("NonLinear converter model is not implemented yet")
    elseif simulation_settings.converter_model == :NF_lossless || simulation_settings.converter_model == :NF_lossy

        if simulation_settings.converter_modularization == :discrete
            
            JuMP.@constraint(model, converter_flow_cap_ac[g in prerequisites_data.conv_ac_side_virtual_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
               grid.Generators[g].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.ac_gen_to_conv_id_dict[g], k] ≤ model[:p_conv_ac][g, k, t] ≤ grid.Generators[g].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.ac_gen_to_conv_id_dict[g], k])

            JuMP.@constraint(model, converter_flow_cap_dc[g in prerequisites_data.conv_dc_side_virtual_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                grid.DCGenerators[g].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.dc_gen_to_conv_id_dict[g], k] ≤ model[:p_conv_dc][g, k, t] ≤ grid.DCGenerators[g].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.dc_gen_to_conv_id_dict[g], k])

        elseif simulation_settings.converter_modularization == :continuous
            Converters = grid.Converters
            
            JuMP.@constraint(model, converter_flow_cap_ac_side_1_up[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:γ_conv][d,k,t]*grid.Generators[Converters[prerequisites_data.conv_duplets[d][1]].gen_ac_id].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][1], k]
                ≤ model[:p_conv_ac][Converters[prerequisites_data.conv_duplets[d][1]].gen_ac_id, k, t])

            JuMP.@constraint(model, converter_flow_cap_ac_side_1_down[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_conv_ac][Converters[prerequisites_data.conv_duplets[d][1]].gen_ac_id, k, t]
                ≤ model[:γ_conv][d,k,t]*grid.Generators[Converters[prerequisites_data.conv_duplets[d][1]].gen_ac_id].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][1], k])
               
            JuMP.@constraint(model, converter_flow_cap_ac_side_2_up[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                (1-model[:γ_conv][d,k,t])*grid.Generators[Converters[prerequisites_data.conv_duplets[d][2]].gen_ac_id].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][2], k]
                ≤ model[:p_conv_ac][Converters[prerequisites_data.conv_duplets[d][2]].gen_ac_id, k, t])

            JuMP.@constraint(model, converter_flow_cap_ac_side_2_down[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_conv_ac][Converters[prerequisites_data.conv_duplets[d][2]].gen_ac_id, k, t] 
                ≤ (1-model[:γ_conv][d,k,t])*grid.Generators[Converters[prerequisites_data.conv_duplets[d][2]].gen_ac_id].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][2], k])
            
            JuMP.@constraint(model, converter_flow_cap_dc_side_1_up[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:γ_conv][d,k,t]*grid.DCGenerators[Converters[prerequisites_data.conv_duplets[d][1]].gen_dc_id].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][1], k]
                ≤ model[:p_conv_dc][Converters[prerequisites_data.conv_duplets[d][1]].gen_dc_id, k, t])

            JuMP.@constraint(model, converter_flow_cap_dc_side_1_down[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_conv_dc][Converters[prerequisites_data.conv_duplets[d][1]].gen_dc_id, k, t]
                ≤ model[:γ_conv][d,k,t]*grid.DCGenerators[Converters[prerequisites_data.conv_duplets[d][1]].gen_dc_id].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][1], k])

            JuMP.@constraint(model, converter_flow_cap_dc_side_2_up[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                (1-model[:γ_conv][d,k,t])*grid.DCGenerators[Converters[prerequisites_data.conv_duplets[d][2]].gen_dc_id].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][2], k]
                ≤ model[:p_conv_dc][Converters[prerequisites_data.conv_duplets[d][2]].gen_dc_id, k, t])

            JuMP.@constraint(model, converter_flow_cap_dc_side_2_down[d in sort(collect(keys(prerequisites_data.conv_duplets))), k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_conv_dc][Converters[prerequisites_data.conv_duplets[d][2]].gen_dc_id, k, t]
                ≤ (1-model[:γ_conv][d,k,t])*grid.DCGenerators[Converters[prerequisites_data.conv_duplets[d][2]].gen_dc_id].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.conv_duplets[d][2], k])
        else
            error("Invalid modularization technique $(modularization)")
        end

        JuMP.@constraint(model, b2b_converter_flow_cap[g in prerequisites_data.b2b_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            grid.Generators[g].Pg_min*prerequisites_data.Contingency_Map["conv"][prerequisites_data.ac_gen_to_conv_id_dict[g], k] ≤ model[:p_conv_b2b][g, k, t] ≤ grid.Generators[g].Pg_max*prerequisites_data.Contingency_Map["conv"][prerequisites_data.ac_gen_to_conv_id_dict[g], k])
    
    # elseif simulation_settings.converter_model == :NF_lossy
    #     error("Linear lossy converter model is not implemented yet")
    #     # JuMP.@constraint(model, conv_power[c in converter_set], p_conv[grid.Converters[c].gen_ac_id] + p_conv[grid.Converters[c].gen_dc_id] == losses)
    else
        error("Requested converter model ($(simulation_settings.converter_model)) is not implemented.")
    end
end

function DOPF_nodal_balance_dc_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.dc_grid_model == :FNL
        error("NonLinear nodal balance is not implemented yet for DC grids")
    elseif simulation_settings.dc_grid_model == :NF
        JuMP.@constraint(model, dc_active_nodal_balance[n in prerequisites_data.dc_node_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            sum([model[:p_gen_dc][g,k,t] for g in grid.DCBuses[n].ConnectedGensIDs if g in prerequisites_data.dc_gen_ids], init=0)
            + sum([model[:p_conv_dc][g,k,t] for g in grid.DCBuses[n].ConnectedGensIDs if g in prerequisites_data.conv_dc_side_virtual_gen_ids], init = 0)
            - sum([grid.DCLoads[d].Pd_t[t] for d in grid.DCBuses[n].ConnectedLoadsIDs if d in prerequisites_data.dc_load_ids], init = 0) 
            == sum([model[:p_branch_dc][l,1,k,t] for l in grid.DCBuses[n].ConnectedLinesIDs if grid.DCBranches[l].Fr_bus_ID == n], init = 0)
            + sum([model[:p_branch_dc][l,2,k,t] for l in grid.DCBuses[n].ConnectedLinesIDs if grid.DCBranches[l].To_bus_ID == n], init = 0)
            )
    else
        error("Requested formulation ($formulation) is not implemented.")
    end
end

function DOPF_powerflow_dc_branch_static!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.dc_grid_model == :NF
        JuMP.@constraint(model, dc_active_flow[l in prerequisites_data.dc_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_branch_dc][l,1,k,t] + model[:p_branch_dc][l,2,k,t] == 0)
    elseif simulation_settings.dc_grid_model == :FNL
        error("Full non-linear model of DC grid is not implemented yet")
    else
        error("$(simulation_settings.dc_grid_model) model of DC grid is not implemented")
    end
end

function DOPF_thermal_limits_dc_branch_static!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.dc_grid_model == :NF
        JuMP.@constraint(model, dc_active_flow_limits[l in prerequisites_data.dc_branch_ids, s in [1,2], k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            -grid.DCBranches[l].rating*prerequisites_data.Contingency_Map["dc_branch"][l,k]*prerequisites_data.base_MVA ≤ model[:p_branch_dc][l,s,k,t] ≤ grid.DCBranches[l].rating*prerequisites_data.Contingency_Map["dc_branch"][l,k]*prerequisites_data.base_MVA)
    elseif simulation_settings.dc_grid_model == :FNL
        error("Full non-linear model of DC grid is not implemented yet")
    else
        error("$(simulation_settings.dc_grid_model) model of DC grid is not implemented")
    end
end

function DOPF_powerflow_dc_link!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.dc_grid_model == :NF
        JuMP.@constraint(model, dc_link_active_flow[l in prerequisites_data.dc_link_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:p_dclink][grid.DCLinks[l].Fr_gen_ID,k,t] + model[:p_dclink][grid.DCLinks[l].To_gen_ID,k,t] == 0)
    elseif simulation_settings.dc_grid_model == :FNL
        error("Full non-linear model of DC grid is not implemented yet")
    else
        error("$(simulation_settings.dc_grid_model) model of DC grid is not implemented")
    end
end

function DOPF_thermal_limits_dc_link!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if simulation_settings.dc_grid_model == :NF
        JuMP.@constraint(model, dc_link_active_flow_limits_fr_side[g in prerequisites_data.dc_link_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            grid.Generators[grid.DCLinks[l].Fr_gen_ID].Pg_min*prerequisites_data.Contingency_Map["dc_link"][l,k] ≤ model[:p_dclink][grid.DCLinks[l].Fr_gen_ID,k,t] ≤ grid.Generators[grid.DCLinks[l].Fr_gen_ID].Pg_max*prerequisites_data.Contingency_Map["dc_link"][l,k])
            
        JuMP.@constraint(model, dc_link_active_flow_limits_to_side[g in prerequisites_data.dc_link_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            grid.Generators[grid.DCLinks[l].To_gen_ID].Pg_min*prerequisites_data.Contingency_Map["dc_link"][l,k] ≤ model[:p_dclink][grid.DCLinks[l].To_gen_ID,k,t] ≤ grid.Generators[grid.DCLinks[l].To_gen_ID].Pg_max*prerequisites_data.Contingency_Map["dc_link"][l,k])
    elseif simulation_settings.dc_grid_model == :FNL
        error("Full non-linear model of DC grid is not implemented yet")
    else
        error("$(simulation_settings.dc_grid_model) model of DC grid is not implemented")
    end
end

function DOPF_voltage_limits_dc_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
end

function DOPF_generator_limits_dc_node!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if length(prerequisites_data.dc_gen_ids) != 0
        JuMP.@constraint(model, dc_gen_cap_non_commitable[g in prerequisites_data.dc_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            grid.DCGenerators[g].P_max*prerequisites_data.Contingency_Map["dc_gen"][g,k] 
            ≤ model[:p_gen_dc][g,k,t] ≤ grid.DCGenerators[g].P_min*prerequisites_data.Contingency_Map["dc_gen"][g,k])
    end
end

function DOPF_transition_constraints!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    ac_gen_id_to_gen_root = prerequisites_data.ac_gen_id_to_gen_root
    # maximum pre-contingency transmission switching
    if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
        if simulation_settings.max_transmission_switching["pre_contingency"] !== Inf
            JuMP.@constraint(model, pre_contingency_maximum_off_ac_branches[t in prerequisites_data.time_horizon],
                sum(1-model[:z_l][l,1,t] for l in prerequisites_data.ac_active_dynamic_branch_ids) ≤ simulation_settings.max_transmission_switching["pre_contingency"])
        end
    end
    # maximum pre-contingency busbar splitting
    if length(prerequisites_data.ac_active_coupler_ids) != 0
        if simulation_settings.max_busbar_splitting["pre_contingency"] !== Inf
            JuMP.@constraint(model, pre_contingency_maximum_off_ac_couplers[t in prerequisites_data.time_horizon],
                sum(1-model[:z_c][l,1,t] for l in prerequisites_data.ac_active_coupler_ids) ≤ simulation_settings.max_busbar_splitting["pre_contingency"])
        end
    end

    # Load shedding allowance
    if simulation_settings.load_shedding == [:post]
        JuMP.@constraint(model, no_load_shedding_pre[d in prerequisites_data.ac_load_shedding_ids, t in prerequisites_data.time_horizon],
            model[:p_ls_ac][d,1,t] == 0)
    elseif simulation_settings.load_shedding == [:pre]
        JuMP.@constraint(model, no_load_shedding_post[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k in prerequisites_data.k_t[t]],
            model[:p_ls_ac][d,k,t] == 0)
    elseif simulation_settings.load_shedding == []
        JuMP.@constraint(model, no_load_shedding_at_all[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k in prerequisites_data.k_t[t]],
            model[:p_ls_ac][d,k,t] == 0)
    end

    # Temporal transitions
    if length(prerequisites_data.time_horizon) ≥ 2
        if length(prerequisites_data.commitable_gen_ids)!=0
            # 1. Unit Commitment logic
            JuMP.@constraint(model, logic[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon[2:end]], model[:u_gt][g,t]-model[:u_gt][g,t-1] == model[:α_gt][g,t] - model[:β_gt][g,t])
            JuMP.@constraint(model, logic_init[g in prerequisites_data.commitable_gen_ids], model[:u_gt][g,1] == model[:α_gt][g,1] + model[:β_gt][g,1])
            
            JuMP.@constraint(model, MUT[g in prerequisites_data.commitable_gen_ids, t in grid.Generators[g].min_up_time:prerequisites_data.time_horizon[end]],
                sum(model[:α_gt][g,i] for i in t-grid.Generators[g].min_up_time+1:t) ≤ model[:u_gt][g,t])
            JuMP.@constraint(model, MDT[g in prerequisites_data.commitable_gen_ids, t in grid.Generators[g].min_down_time:prerequisites_data.time_horizon[end]],
                sum(model[:β_gt][g,i] for i in t-grid.Generators[g].min_down_time+1:t) ≤ 1-model[:u_gt][g,t])
                
            # 2.1 Temporal Scheduling
            JuMP.@constraint(model, Ramp_Up_commitable[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t]+prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t]-model[:p_gen_ac][g,1,t-1]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t-1] ≤ grid.Generators[g].Δ_up)
            JuMP.@constraint(model, Ramp_Down_commitable[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t-1]+prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t-1]-model[:p_gen_ac][g,1,t]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t] ≤ grid.Generators[g].Δ_down)
        end

        if length(prerequisites_data.non_commitable_gen_ids) != 0
            # 2.2 Temporal Scheduling
            JuMP.@constraint(model, Ramp_Up_non_commitable[g in prerequisites_data.non_commitable_gen_ids, t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t]+prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t]-model[:p_gen_ac][g,1,t-1]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t-1] ≤ grid.Generators[g].Δ_up)
            JuMP.@constraint(model, Ramp_Down_non_commitable[g in prerequisites_data.non_commitable_gen_ids, t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t-1]+prerequisites_data.Schedule["gen"][g][t-1]-model[:p_gen_ac][ac_gen_id_to_gen_root[g],1,t]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t] ≤ grid.Generators[g].Δ_down)
        end

        if length(keys(prerequisites_data.fixed_commitments)) != 0
            # 2.3 Temporal Scheduling
            JuMP.@constraint(model, Ramp_Up_fixed_commitments[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t]+prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t]-model[:p_gen_ac][g,1,t-1]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t-1] ≤ grid.Generators[g].Δ_up)
            JuMP.@constraint(model, Ramp_Down_fixed_commitments[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon[2:end]], model[:p_gen_ac][g,1,t-1]+prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t-1]-model[:p_gen_ac][g,1,t]-prerequisites_data.Schedule["gen"][ac_gen_id_to_gen_root[g]][t] ≤ grid.Generators[g].Δ_down)
        end

        if simulation_settings.activate_temporal_switching
            if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
                # 3.1 Temporal Transmission Switching --> minimum cool-down time (after changing status), maximum switching per step
                MCDT_l = simulation_settings.max_transmission_switching["MCDT"]
                JuMP.@constraint(model, logic_l[l in prerequisites_data.ac_active_dynamic_branch_ids, t in prerequisites_data.time_horizon[2:end]], model[:z_l][l,1,t]-model[:z_l][l,1,t-1] == model[:α_l][l,t] - model[:β_l][l,t])
                JuMP.@constraint(model, logic_init_l[l in prerequisites_data.ac_active_dynamic_branch_ids], model[:z_l][l,1,1] == model[:α_l][l,1] + model[:β_l][l,1])
                
                JuMP.@constraint(model, MUT_l[l in prerequisites_data.ac_active_dynamic_branch_ids, t in MCDT_l:prerequisites_data.time_horizon[end]],
                    sum(model[:α_l][l,i] for i in t-MCDT_l+1:t) ≤ model[:z_l][l,1,t])
                JuMP.@constraint(model, MDT_l[l in prerequisites_data.ac_active_dynamic_branch_ids, t in MCDT_l:prerequisites_data.time_horizon[end]],
                    sum(model[:β_l][l,i] for i in t-MCDT_l+1:t) ≤ 1-model[:z_l][l,1,t])
            end

            if length(prerequisites_data.ac_active_reconf_ids) != 0
                # 3.2 Temporal Substation Reconfiguration --> minimum cool-down time (after changing status)
                MCDT_r = simulation_settings.max_substation_reconf["MCDT"]
                JuMP.@constraint(model, logic_r[l in prerequisites_data.ac_active_reconf_ids, t in prerequisites_data.time_horizon[2:end]], model[:z_r][l,1,t]-model[:z_r][l,1,t-1] == model[:α_r][l,t] - model[:β_r][l,t])
                JuMP.@constraint(model, logic_init_r[l in prerequisites_data.ac_active_reconf_ids], model[:z_r][l,1,1] == model[:α_r][l,1] + model[:β_r][l,1])
                
                JuMP.@constraint(model, MUT_r[l in prerequisites_data.ac_active_reconf_ids, t in MCDT_r:prerequisites_data.time_horizon[end]],
                    sum(model[:α_r][l,i] for i in t-MCDT_r+1:t) ≤ model[:z_r][l,1,t])
                JuMP.@constraint(model, MDT_r[l in prerequisites_data.ac_active_reconf_ids, t in MCDT_r:prerequisites_data.time_horizon[end]],
                    sum(model[:β_r][l,i] for i in t-MCDT_r+1:t) ≤ 1-model[:z_r][l,1,t])
            end

            if length(prerequisites_data.ac_active_coupler_ids) != 0
                # 3.3 Temporal Busbar Splitting --> minimum cool-down time (after changing status), maximum splitting per step
                MCDT_c = simulation_settings.max_substation_reconf["MCDT"]
                JuMP.@constraint(model, logic_c[l in prerequisites_data.ac_active_coupler_ids, t in prerequisites_data.time_horizon[2:end]], model[:z_c][l,1,t]-model[:z_c][l,1,t-1] == model[:α_c][l,t] - model[:β_c][l,t])
                JuMP.@constraint(model, logic_init_c[l in prerequisites_data.ac_active_coupler_ids], model[:z_c][l,1,1] == model[:α_c][l,1] + model[:β_c][l,1])
                
                JuMP.@constraint(model, MUT_c[l in prerequisites_data.ac_active_coupler_ids, t in MCDT_c:prerequisites_data.time_horizon[end]],
                    sum(model[:α_c][l,i] for i in t-MCDT_c+1:t) ≤ model[:z_c][l,1,t])
                JuMP.@constraint(model, MDT_c[l in prerequisites_data.ac_active_coupler_ids, t in MCDT_c:prerequisites_data.time_horizon[end]],
                    sum(model[:β_c][l,i] for i in t-MCDT_c+1:t) ≤ 1-model[:z_c][l,1,t])
            end
        end
    end

    # Contingency transitions
    if length(prerequisites_data.k) ≥ 2

        # 1.1 maximum post-contingency transmission switching
        if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
            if simulation_settings.max_transmission_switching["post_contingency"] !== Inf
                JuMP.@constraint(model, post_contingency_maximum_off_ac_branches[k in prerequisites_data.k[2:end],t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                    sum(1-model[:z_l][l,k,t] for l in prerequisites_data.ac_active_dynamic_branch_ids) ≤ simulation_settings.max_transmission_switching["post_contingency"])
            end
        end
        # 1.2 maximum post-contingency busbar splitting
        if length(prerequisites_data.ac_active_coupler_ids) != 0
            if simulation_settings.max_busbar_splitting["post_contingency"] !== Inf       
                JuMP.@constraint(model, post_contingency_maximum_off_ac_couplers[k in prerequisites_data.k[2:end],t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                    sum(1-model[:z_c][l,k,t] for l in prerequisites_data.ac_active_coupler_ids) ≤ simulation_settings.max_busbar_splitting["post_contingency"])
            end
        end

        # 2 Corrective redispatch constraints
        if length(prerequisites_data.contingency_redispatch) == 0
            JuMP.@constraint(model,no_redispatch[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,1,t] == model[:p_gen_ac][g,k,t])
        else
            JuMP.@constraint(model,Ramp_Up_contingency_redispatch[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.contingency_redispatch, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                (model[:p_gen_ac][g,1,t] - model[:p_gen_ac][g,k,t])*prerequisites_data.Contingency_Map["ac_gen"][g,k] ≤ grid.Generators[g].Δ_down)
            JuMP.@constraint(model,Ramp_Down_contingency_redispatch[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.contingency_redispatch, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]], 
                (model[:p_gen_ac][g,k,t] - model[:p_gen_ac][g,1,t])*prerequisites_data.Contingency_Map["ac_gen"][g,k] ≤ grid.Generators[g].Δ_up)
                
            fixed_redispatch_k = collect(setdiff(setdiff(Set(prerequisites_data.k),Set(prerequisites_data.contingency_redispatch)), Set([1])))
            JuMP.@constraint(model,no_redispatch[g in prerequisites_data.ac_gen_ids, k in fixed_redispatch_k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:p_gen_ac][g,1,t] == model[:p_gen_ac][g,k,t])
        end

        # 3 Dynamic converter control constraints
        if ! simulation_settings.dynamic_converter_control
            if length(prerequisites_data.conv_ids) != 0
                JuMP.@constraint(model,static_converters_ac[g in prerequisites_data.conv_ac_side_virtual_gen_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                    model[:p_conv_ac][g,1,t] == model[:p_conv_ac][g,k,t])
                JuMP.@constraint(model,static_converters_dc[g in prerequisites_data.conv_dc_side_virtual_gen_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                    model[:p_conv_dc][g,1,t] == model[:p_conv_dc][g,k,t])
            end
            if length(prerequisites_data.dc_link_ids) != 0
                JuMP.@constraint(model,static_dc_link[g in prerequisites_data.b2b_gen_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                    model[:p_conv_b2b][g,1,t] == model[:p_conv_b2b][g,k,t])
            end
        end

        # 4 Corrective transmission switching constraints
        if simulation_settings.transmission_switching == [:post]
            JuMP.@constraint(model, static_grid_pre_contingency[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:z_l][l,1,t] == 1)
        elseif simulation_settings.transmission_switching == [:pre]
            JuMP.@constraint(model, static_grid_post_contingency[l in prerequisites_data.ac_active_dynamic_branch_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:z_l][l,k,t] == model[:z_l][l,1,t])
        end

        # 5 Corrective substation switching constraints
        if simulation_settings.substation_switching["reconf"] == [:pre]
            JuMP.@constraint(model, static_substation_post_contingency[r in prerequisites_data.ac_active_reconf_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:z_r][r,k,t] == model[:z_r][r,1,t])
        elseif simulation_settings.substation_switching["reconf"] == [:post]
            JuMP.@constraint(model, static_substation_pre_contingency[t in prerequisites_data.time_horizon], sum([model[:z_r][l,1,t] for l in prerequisites_data.normally_opened_reconf_lines], init=0) == 0)
        end

        if simulation_settings.substation_switching["splitting"] == [:pre]
            JuMP.@constraint(model, static_coupler_post_contingency[c in prerequisites_data.ac_active_coupler_ids, k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
                model[:z_c][c,k,t] == model[:z_c][c,1,t])
        elseif simulation_settings.substation_switching["splitting"] == [:post]
            JuMP.@constraint(model, static_coupler_pre_contingency[c in prerequisites_data.ac_active_coupler_ids, t in prerequisites_data.time_horizon], model[:z_c][c,1,t] == 1)
        end

        # 6 Substation switching novel bound-tightening constraints
        if simulation_settings.substation_switching["reconf"] == [:pre] && simulation_settings.substation_switching["splitting"] == [:pre]
            JuMP.@constraint(model,switched_off_reconf_lines_pre[s in keys(grid.Substations), t in prerequisites_data.time_horizon; length(grid.Substations[s].Reconf_CouplerLines_IDs) == 1], 
                sum(model[:z_r][l,1,t] for l in grid.Substations[s].Reconf_AuxLines_IDs if l in prerequisites_data.normally_opened_reconf_lines) ≤ (1-model[:z_c][grid.Substations[s].Reconf_CouplerLines_IDs[1], 1, t])*length(grid.Substations[s].Reconf_AuxLines_IDs)/2)
        elseif simulation_settings.substation_switching["reconf"] == [:post] && simulation_settings.substation_switching["splitting"] == [:post]
            JuMP.@constraint(model,switched_off_reconf_lines_post[s in keys(grid.Substations), k in prerequisites_data.k[2:end], t in prerequisites_data.time_horizon; length(grid.Substations[s].Reconf_CouplerLines_IDs) == 1 && k ∈ prerequisites_data.k_t[t]], 
                sum(model[:z_r][l,k,t] for l in grid.Substations[s].Reconf_AuxLines_IDs if l in prerequisites_data.normally_opened_reconf_lines) ≤ (1-model[:z_c][grid.Substations[s].Reconf_CouplerLines_IDs[1], k, t])*length(grid.Substations[s].Reconf_AuxLines_IDs)/2)
        end

    end

end

function DOPF_schedule_fixes!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites; k=1)
    if length(keys(prerequisites_data.fixed_schedules)) != 0
        if simulation_settings.ac_grid_model == :Bθ
            JuMP.@constraint(model, fixed_ac_gen_schedule[g in keys(prerequisites_data.fixed_schedules), t in prerequisites_data.time_horizon], 
                model[:p_gen_ac][g,k,t] == prerequisites_data.fixed_schedules[g][t])
        elseif simulation_settings.ac_grid_model == :AC
        end
    end
end

function DOPF_topology_fixes!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if length(prerequisites_data.ac_fixed_dynamic_branch_ids) != 0
        JuMP.@constraint(model, fixed_transmission_status[l in prerequisites_data.ac_fixed_dynamic_branch_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:z_l_f][l,k,t] == prerequisites_data.fixed_topology[l][k][t])
    end

    if length(prerequisites_data.ac_fixed_reconf_ids) != 0
        JuMP.@constraint(model, fixed_reconf_status[l in prerequisites_data.ac_fixed_reconf_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:z_r_f][l,k,t] == prerequisites_data.fixed_topology[l][k][t])
    end

    if length(prerequisites_data.ac_fixed_coupler_ids) != 0
        JuMP.@constraint(model, fixed_coupler_status[l in prerequisites_data.ac_fixed_coupler_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
            model[:z_c_f][l,k,t] == prerequisites_data.fixed_topology[l][k][t])
    end
end

function DOPF_commitment_fixes!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    if length(keys(prerequisites_data.fixed_commitments)) != 0
        JuMP.@constraint(model, fixed_commitment_status[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon],
            model[:u_gt_f][g,t] == prerequisites_data.fixed_commitments[g][t])
    end
end

function DOPF_objective_function!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    # non_commitable_gen costs + commitable_gen costs + fixed_commitment_gen costs + fixed_schedule_gen costs + load_shedding costs + dc_gen cost
    GenBids = prerequisites_data.Order_Book.Gen_bids
    LoadBids = prerequisites_data.Order_Book.Load_bids
    JuMP.@objective(model, Min, sum([model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in prerequisites_data.non_commitable_gen_ids, t in prerequisites_data.time_horizon], init=0)
        + sum([model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] + grid.Generators[g].C0*model[:u_gt][g,t]+model[:α_gt][g,t]*grid.Generators[g].start_up_cost+model[:β_gt][g,t]*grid.Generators[g].shut_down_cost for g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], init=0)
        + sum([model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon], init=0)
        + sum([model[:p_gen_ac][g,1,t]*GenBids[g]["price"][t][1] for g in keys(prerequisites_data.fixed_schedules), t in prerequisites_data.time_horizon], init=0)
        + sum([model[:p_ls_ac][d,k,t]*LoadBids[d]["price"][t][1] for d  in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon if k ∈ prerequisites_data.k_t[t]], init=0)
        + sum([model[:p_gen_dc][g,1,t]*grid.DCGenerators[g].C1 for g in prerequisites_data.dc_gen_ids, t in prerequisites_data.time_horizon], init=0))
end

###############################################

function DOPF_single_node_balance!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    """
    Don't forget to update `prerequisites_data.k` , `prerequisites_data.k_t`, and `prerequisites_data.time_horizon` before using this if needed!!
    """
    JuMP.@constraint(model, single_node_balance[k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]],
        sum([model[:p_gen_ac][g,k,t] for g in prerequisites_data.ac_gen_ids], init=0)
        + sum([get(get(prerequisites_data.Schedule["gen"],g,Dict()),t,0) for g in keys(prerequisites_data.root_gen_to_duplicate_gen)], init=0)
        + sum([model[:p_gen_dc][g,k,t] for g in prerequisites_data.dc_gen_ids], init=0)
        + sum([model[:p_ls_ac][d,k,t] for d in prerequisites_data.ac_load_shedding_ids], init=0)
        == sum([grid.Loads[d].Pd_t[t] for d in prerequisites_data.ac_load_ids], init = 0)
        + sum([grid.DCLoads[d].Pd_t[t] for d in prerequisites_data.dc_load_ids], init = 0))
end

function DOPF_single_node_variable_initialization!(model::Model, grid ::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    """
    Single node model variable initialization for DOPF simulation, based on simulation settings `simulation_settings` and prerequisites `prerequisites_data`
    """
    # variable initialization
    if simulation_settings.ac_grid_model == :Bθ

        ###################################################################################
        JuMP.@variable(model, p_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])
        JuMP.@variable(model, p_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])

        ###################################################################################
        # Binary commitment variables
        if length(prerequisites_data.commitable_gen_ids) != 0
            JuMP.@variable(model, u_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin)
            JuMP.@variable(model, α_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # start-up
            JuMP.@variable(model, β_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # shut-down
        end

        # Fixed commitment variables
        if length(keys(prerequisites_data.fixed_commitments)) != 0
            JuMP.@variable(model, u_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon])
        end
        ###################################################################################

    elseif simulation_settings.ac_grid_model == :AC
        JuMP.@variable(model, p_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        JuMP.@variable(model, q_gen_ac[g in prerequisites_data.ac_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        
        if length(prerequisites_data.ac_load_shedding_ids) != 0
            JuMP.@variable(model, p_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
            JuMP.@variable(model, q_ls_ac[d in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon])
        end
        ###################################################################################
        # Binary commitment variables
        if length(prerequisites_data.commitable_gen_ids) != 0
            JuMP.@variable(model, u_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin)
            JuMP.@variable(model, α_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # start-up
            JuMP.@variable(model, β_gt[g in prerequisites_data.commitable_gen_ids, t in prerequisites_data.time_horizon], Bin) # shut-down
        end

        # Fixed commitment variables
        if length(keys(prerequisites_data.fixed_commitments)) != 0
            JuMP.@variable(model, u_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon])
            JuMP.@variable(model, α_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon]) # start-up
            JuMP.@variable(model, β_gt_f[g in keys(prerequisites_data.fixed_commitments), t in prerequisites_data.time_horizon]) # shut-down
        end
        
    else
        error("Invalid AC grid model setting ($(simulation_settings.ac_grid_model))")
        return -1
    end

    JuMP.@variable(model, p_gen_dc[g in prerequisites_data.dc_gen_ids, k in prerequisites_data.k, t in prerequisites_data.time_horizon; k ∈ prerequisites_data.k_t[t]])

end