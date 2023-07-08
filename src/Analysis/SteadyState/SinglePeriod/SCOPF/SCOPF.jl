include("SCOPFConstraints.jl")

function build_SCOPF_model!(grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites; optimize=false, prioritize_splitting=false)

    # Sanity check
    solver = simulation_settings.MILP_solver
    NLP_flag = false
    if simulation_settings.ac_grid_model == :ACOPF || simulation_settings.dc_grid_model == :Nonlinear
        solver = simulation_settings.NLP_solver
        NLP_flag = true
    end

    if NLP_flag
        if simulation_settings.transmission_switching != [] || simulation_settings.substation_switching != []
            error("Switching is not supported with Nonlinear problem settings.")
            return -1
        end
    end

    # Model initialization
    model = Model(solver)
    # model = direct_model(solver())

    # Variables initialization
    single_period_SCOPF_variable_initialization!(model, simulation_settings, prerequisites_data)

    # AC grid constraints:
    single_period_SC_angle_limits_ac_grid!(prerequisites_data, grid, model)
    single_period_SC_voltage_limits_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_generator_limits_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_nodal_balance_ac_node!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    if simulation_settings.transmission_switching != []
        single_period_SC_switched_transmission_capacity_limits_ac_branch!(prerequisites_data, grid, model)
        single_period_SC_switched_powerflow_ac_branch!(prerequisites_data, grid, model; max_op=simulation_settings.max_transmission_switching)
    end

    single_period_SC_transmission_capacity_limits_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_powerflow_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    if prerequisites_data.coupler_set != []
        single_period_SC_reconf_split_constraints_ac_grid!(prerequisites_data, grid, model; max_reconf=simulation_settings.max_substation_reconf, max_splitting=simulation_settings.max_busbar_splitting)
    end

    # DC grid constraints:
    if length(keys(grid.DCBuses)) != 0
        single_period_SC_voltage_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_powerflow_dc_branch!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_transmission_capacity_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_nodal_balance_dc_node!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
    end

    # Converter constraints:
    if length(keys(grid.Converters)) != 0 || length(keys(grid.DCLinks)) != 0
        single_period_SC_converter_constraints!(prerequisites_data, grid, model, simulation_settings.converter_model)
        single_period_SC_converter_flow!(prerequisites_data, grid, model, simulation_settings.converter_model; modularization = simulation_settings.converter_modularization)
    end

    single_period_SC_corrective_measures!(prerequisites_data, grid, model, simulation_settings)

    transmission_switching = prerequisites_data.switched_transmission_set != [] ? true : false
    substation_switching = prerequisites_data.coupler_set != [] ? true : false
    single_period_SC_objective!(prerequisites_data, grid, model, transmission_switching, substation_switching)

    if prerequisites_data.coupler_set != [] && prioritize_splitting
        for c in prerequisites_data.coupler_set
            for k in prerequisites_data.k
                MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), model[:z_c][c, k], 1)
            end
        end
    end

    if optimize
        optimize!(model)
        result = SCOPF_post_processing!(grid, model, simulation_settings, prerequisites_data)
        if result == -1
            return -1
        end
        # reduce_grid_SC!(grid,clean=clean)
        update_grid_tables_SC!(grid)
    end

    return model
end

function update_grid_tables_SC!(grid ::PowerGrid;k=1)

    grid.BusData_output = DataFrame(Bus = Int64[], V = Float64[], δ = Float64[],
        Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[],type=[])
    
    grid.LineLoading = DataFrame(BranchID = Int64[], FromBus = Int64[],ToBus = Int64[]
        ,PL_1 = Float64[],PL_2 = Float64[],
        PLoss = Float64[],QL_1 = Float64[]
        ,QL_2 = Float64[],QLoss = Float64[],Utilization = Float64[],type=[])

    grid.Converter_flow = DataFrame(ConverterID = Int64[], AC_Bus = Int64[], DC_Bus = Int64[],
        P_ACDC = Float64[], P_DCAC = Float64[], PLoss = Float64[], Utilization = Float64[],type = [])

    # Populate "BusData_output" dataframe
    # Bus = Int64[], V = Float64[], δ = Float64[],Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[]
    for bus_id in sort(collect(keys(grid.Buses)))
        my_bus = grid.Buses[bus_id]
        condition_1 = my_bus.ConnectedGensIDs != []
        condition_2 = Set([:virtual]) != Set([grid.Generators[g].GenType for g in my_bus.ConnectedGensIDs])
        condition = condition_1 && condition_2
        Pg = (condition) ? sum(grid.Generators[g].Pg_k[k] for g in my_bus.ConnectedGensIDs if grid.Generators[g].GenType != :virtual) : 0
        Qg = (condition) ? sum(grid.Generators[g].Qg_k[k] for g in my_bus.ConnectedGensIDs if grid.Generators[g].GenType != :virtual) : 0
        Pd = (length(my_bus.ConnectedLoadsIDs) != 0) ? sum(grid.Loads[d].Pd for d in my_bus.ConnectedLoadsIDs) : 0
        Qd = (length(my_bus.ConnectedLoadsIDs) != 0) ? sum(grid.Loads[d].Qd for d in my_bus.ConnectedLoadsIDs) : 0
        push!(grid.BusData_output,[bus_id, my_bus.V_magnitude_k[k], my_bus.δ_k[k], Pg, Qg, Pd, Qd, "AC"])
    end

    for bus_id in sort(collect(keys(grid.DCBuses)))
        my_bus = grid.DCBuses[bus_id]
        push!(grid.BusData_output,[bus_id, my_bus.V_magnitude_k[k], 0, 0, 0, 0, 0, "DC"])
    end

    # Populate "LineLoading" dataframe
    # BranchID = Int64[], FromBus = Int64[],ToBus = Int64[],PL_1 = Float64[],PL_2 = Float64[],PLoss = Float64[],QL_1 = Float64[],QL_2 = Float64[],QLoss = Float64[],Utilization = Float64[]
    for line_id in sort(collect(keys(grid.Branches)))
        my_branch = grid.Branches[line_id]
        FromBus = my_branch.Fr_bus_ID
        ToBus = my_branch.To_bus_ID
        PL_1 = my_branch.PowerFlow_ij_k[k]
        PL_2 = my_branch.PowerFlow_ji_k[k]
        PLoss = my_branch.losses_P_k[k]
        QL_1 = my_branch.ReactFlow_ij_k[k]
        QL_2 = my_branch.ReactFlow_ji_k[k]
        QLoss = my_branch.losses_Q_k[k]
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating*grid.S_base), sqrt(PL_2^2+QL_2^2)/(my_branch.rating*grid.S_base)])
        if my_branch.BranchType == 0
            label = "AC"
        elseif my_branch.BranchType == 1
            label = "RECONF"
        elseif my_branch.BranchType == 2
            label = "COUPLER"
        end
        push!(grid.LineLoading,[line_id,FromBus,ToBus,PL_1,PL_2,PLoss,QL_1,QL_2,QLoss,Utilization*100, label])
    end

    for line_id in sort(collect(keys(grid.DCBranches)))
        my_branch = grid.DCBranches[line_id]
        FromBus = my_branch.Fr_bus_ID
        ToBus = my_branch.To_bus_ID
        PL_1 = my_branch.PowerFlow_ij_k[k]
        PL_2 = my_branch.PowerFlow_ji_k[k]
        PLoss = my_branch.losses_P_k[k]
        QL_1 = 0
        QL_2 = 0
        QLoss = 0
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating*grid.S_base), sqrt(PL_2^2+QL_2^2)/(my_branch.rating*grid.S_base)])
        push!(grid.LineLoading,[line_id,FromBus,ToBus,PL_1,PL_2,PLoss,QL_1,QL_2,QLoss,Utilization*100, "DC"])
    end

    for line_id in sort(collect(keys(grid.DCLinks)))
        my_branch = grid.DCLinks[line_id]
        FromBus = my_branch.Fr_bus_ID
        ToBus = my_branch.To_bus_ID
        PL_1 = my_branch.PowerFlow_ij_k[k]
        PL_2 = my_branch.PowerFlow_ji_k[k]
        PLoss = my_branch.losses_P_k[k]
        QL_1 = 0
        QL_2 = 0
        QLoss = 0
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating), sqrt(PL_2^2+QL_2^2)/(my_branch.rating)])
        push!(grid.LineLoading,[line_id,FromBus,ToBus,PL_1,PL_2,PLoss,QL_1,QL_2,QLoss,Utilization*100, "DCLink"])
    end

    # Flow through Converters
    for conv in sort(collect(keys(grid.Converters)))
        ConverterID = grid.Converters[conv].Conv_ID
        AC_Bus = grid.Generators[grid.Converters[conv].gen_ac_id].GenBus_ID
        DC_Bus = grid.Generators[grid.Converters[conv].gen_dc_id].GenBus_ID
        P_ACDC = grid.Generators[grid.Converters[conv].gen_ac_id].Pg_k[k]
        P_DCAC = grid.Generators[grid.Converters[conv].gen_dc_id].Pg_k[k]
        PLoss = abs(P_ACDC + P_DCAC)
        Utilization = maximum([abs(P_ACDC)/grid.Converters[conv].rate,abs(P_ACDC)/grid.Converters[conv].rate])
        type = string(grid.Converters[conv].type)
        push!(grid.Converter_flow, [ConverterID, AC_Bus, DC_Bus, P_ACDC,P_DCAC,PLoss,Utilization*100,type])
    end
end

function compile_simulation_prerequisites_SC!(grid ::PowerGrid, contingency_type ::Symbol,include_leafs, include_HVDC,
         switchable_lines ::Array, substations ::Array,
         hybrid_substations ::Array,B2B_capacities ::Dict;reference_node=nothing)
    
    activate_transmission_switch!(grid , switchable_lines)
    convert_bus2substation!(grid,substations,false) # Regular substations
    convert_bus2substation!(grid,hybrid_substations,true,B2B_capacities) # Soft switch substations

    ac_Nodes_set = collect(Int.(keys(grid.Buses)))
    aux_bus_set = [key for key in ac_Nodes_set if grid.Buses[key].BusType == 1]
    dc_Nodes_set = collect(Int.(keys(grid.DCBuses)))

    Branch_set = collect(Int.(keys(grid.Branches)))
    Branch_nodes = [Set([grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID]) for branch in Branch_set]
    branch_dictionary = Dict()
    for branch_id in Branch_set
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    Transmission_set = [key for key in Branch_set if grid.Branches[key].BranchType == 0]
    dc_Transmission_set = keys(grid.DCBranches)
    unswitched_Transmission_set = collect(setdiff(Set(Transmission_set),Set(switchable_lines)))
    unswitched_Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in unswitched_Transmission_set]
    switched_Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in switchable_lines]
    
    dc_branch_dicitionary = Dict()
    for branch_id in dc_Transmission_set
        push!(dc_branch_dicitionary, Set([grid.DCBranches[branch_id].Fr_bus_ID,grid.DCBranches[branch_id].To_bus_ID]) => branch_id)
    end

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

    Transmission_nodes = [Set([grid.Branches[key].Fr_bus_ID,grid.Branches[key].To_bus_ID]) for key in Transmission_set]
    dc_Transmission_nodes = [Set([grid.DCBranches[key].Fr_bus_ID,grid.DCBranches[key].To_bus_ID]) for key in dc_Transmission_set]
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
    load_set = keys(grid.Loads)

    B = -1*imag(grid.Y_bus)

    ############################### Converters ###############################
    ac_virtual_gen_set = []
    for conv_id in keys(grid.Converters)
        if grid.Converters[conv_id].type == :ACDC
            push!(ac_virtual_gen_set,grid.Converters[conv_id].gen_ac_id)
        elseif grid.Converters[conv_id].type == :B2B
            push!(ac_virtual_gen_set,grid.Converters[conv_id].gen_ac_id)
            push!(ac_virtual_gen_set,grid.Converters[conv_id].gen_dc_id)
        end
    end

    dc_virtual_gen_set = []
    for conv_id in keys(grid.Converters)
        if grid.Converters[conv_id].type == :ACDC
            push!(dc_virtual_gen_set,grid.Converters[conv_id].gen_dc_id)
        end
    end
    ########################################################################

    ############################### DC Link ###############################
    dc_link_virtual_gen_set = []
    for link_id in keys(grid.DCLinks)
        push!(dc_link_virtual_gen_set,grid.DCLinks[link_id].Fr_gen_ID)
        push!(dc_link_virtual_gen_set,grid.DCLinks[link_id].To_gen_ID)
    end
    ########################################################################
    converter_set = keys(grid.Converters)
    default_off_reconf = [key for key in Reconf_set if grid.Branches[key].GeneralSwitch.SwitchingStatus == 0]
    all_gen_set = keys(grid.Generators)

    b2b_gen_set = []
    b2b_coupler_dict = Dict()
    for conv_id in keys(grid.Converters)
        if grid.Converters[conv_id].type == :B2B
            push!(b2b_gen_set,grid.Converters[conv_id].gen_dc_id)
            push!(b2b_gen_set,grid.Converters[conv_id].gen_ac_id)
            coupler_id = Coupler_dict[Set([grid.Converters[conv_id].DC_Bus_ID,grid.Converters[conv_id].AC_Bus_ID])]
            push!(b2b_coupler_dict,grid.Converters[conv_id].gen_dc_id => coupler_id)
            push!(b2b_coupler_dict,grid.Converters[conv_id].gen_ac_id => coupler_id)
        end
    end

    if contingency_type == :Transmission_Contingencies
        a_g, a_l, a_l_dc_branch, a_l_dc_link, k = Generate_Contingencies!(grid,:TC; include_leafs=include_leafs, include_HVDC=include_HVDC)
    elseif contingency_type == :TransGen_Contingencies
        a_g, a_l, a_l_dc_branch, a_l_dc_link, k = Generate_Contingencies!(grid,:TGC; include_leafs=include_leafs, include_HVDC=include_HVDC)
    else 
        error("Invalid contingency type ($contingency_type)")
    end

    k = 1:k

    scopf_prereq = SCOPF_Prerequisites(ac_Nodes_set,aux_bus_set,Branch_nodes,Gen_set,load_set,
        ac_virtual_gen_set,dc_virtual_gen_set,all_gen_set,B,Sbase,unswitched_Transmission_nodes,branch_dictionary,
        switched_Transmission_nodes,Reconf_nodes,Reconf_dict,Coupler_nodes,Coupler_dict,default_off_reconf,Coupler_set,
        converter_set,dc_Nodes_set,dc_Transmission_nodes,dc_branch_dicitionary,b2b_gen_set,b2b_coupler_dict,switchable_lines,Reconf_set,a_g, a_l, a_l_dc_branch, a_l_dc_link, k,reference_node)

    return scopf_prereq
end

function SCOPF_post_processing!(grid ::PowerGrid, model ::Model, simulation_settings ::SCOPF_SimulationSettings, prerequisites ::SCOPF_Prerequisites)
    
    if ! JuMP.has_values(model)
        return -1
    end

    grid.Operating_Cost = JuMP.objective_value(model)
    Branch_set = keys(grid.Branches)
    dc_Branch_set = keys(grid.DCBranches)
    k = sort(collect(prerequisites.k))
    # AC grid data
    if simulation_settings.ac_grid_model == :DCOPF
        δ = JuMP.value.(model[:δ])
        pij = JuMP.value.(model[:pij])
        p = JuMP.value.(model[:p])
        
        [grid.Generators[g].Pg = p[g,1] for g in prerequisites.gen_set]
        [grid.Generators[g].Qg = 0 for g in prerequisites.gen_set]
        [grid.Buses[bus].V_magnitude = 1 for bus in prerequisites.nodes_set]
        [grid.Buses[bus].δ = δ[bus,1] for bus in prerequisites.nodes_set]

        [grid.Generators[g].Pg_k = [p[g,κ] for κ in k] for g in prerequisites.gen_set]
        [grid.Generators[g].Qg_k = zeros(length(prerequisites.k)) for g in prerequisites.gen_set]
        [grid.Buses[bus].V_magnitude_k = ones(length(prerequisites.k)) for bus in prerequisites.nodes_set]
        [grid.Buses[bus].δ_k = [δ[bus,κ] for κ in k] for bus in prerequisites.nodes_set]

        [grid.Branches[branch].PowerFlow_ij_k = [pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji_k = [pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].losses_P_k = abs.(grid.Branches[branch].PowerFlow_ij_k+grid.Branches[branch].PowerFlow_ji_k) for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij_k = zeros(length(prerequisites.k)) for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji_k = zeros(length(prerequisites.k)) for branch in Branch_set]
        [grid.Branches[branch].losses_Q_k = zeros(length(prerequisites.k)) for branch in Branch_set]
        

    elseif simulation_settings.ac_grid_model == :ACOPF
        
        v = JuMP.value.(model[:v])
        δ = JuMP.value.(model[:δ])
        pij = JuMP.value.(model[:pij])
        qij = JuMP.value.(model[:qij])
        p = JuMP.value.(model[:p])
        q = JuMP.value.(model[:q])

        [grid.Generators[g].Pg_k = [p[g,κ] for κ in k] for g in prerequisites.gen_set]
        [grid.Generators[g].Qg_k = [q[g,κ] for κ in k] for g in prerequisites.gen_set]
        [grid.Buses[bus].V_magnitude_k = [v[bus,κ] for κ in k] for bus in prerequisites.nodes_set]
        [grid.Buses[bus].δ_k = [δ[bus,κ] for κ in k ]  for bus in prerequisites.nodes_set]

        [grid.Branches[branch].PowerFlow_ij_k = [pij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].PowerFlow_ji_k = [pij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].losses_P_k = abs.(grid.Branches[branch].PowerFlow_ij_k+grid.Branches[branch].PowerFlow_ji_k) for branch in Branch_set]

        [grid.Branches[branch].ReactFlow_ij_k = [qij[grid.Branches[branch].Fr_bus_ID,grid.Branches[branch].To_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].ReactFlow_ji_k = [qij[grid.Branches[branch].To_bus_ID,grid.Branches[branch].Fr_bus_ID,κ] for κ in k] for branch in Branch_set]
        [grid.Branches[branch].losses_Q_k = abs.(grid.Branches[branch].ReactFlow_ij_k+grid.Branches[branch].ReactFlow_ji_k) for branch in Branch_set]
    else
        error("Unidentified grid model setting ($simulation_settings.ac_grid_model)")
        return -1
    end

    # DC grid data
    if simulation_settings.dc_grid_model == :Linear
        if length(collect(dc_Branch_set)) != 0
            pij_dc = JuMP.value.(model[:pij_dc])
            [grid.DCBranches[branch].PowerFlow_ij_k = [pij_dc[grid.DCBranches[branch].Fr_bus_ID,grid.DCBranches[branch].To_bus_ID,κ] for κ in k] for branch in dc_Branch_set]
            [grid.DCBranches[branch].PowerFlow_ji_k = [pij_dc[grid.DCBranches[branch].To_bus_ID,grid.DCBranches[branch].Fr_bus_ID,κ] for κ in k] for branch in dc_Branch_set]
            [grid.DCBranches[branch].losses_P_k = abs.(grid.DCBranches[branch].PowerFlow_ij_k+grid.DCBranches[branch].PowerFlow_ji_k) for branch in dc_Branch_set]
        end
        if length(prerequisites.dc_Nodes_set) != 0
            [grid.DCBuses[bus].V_magnitude_k = ones(length(prerequisites.k)) for bus in prerequisites.dc_Nodes_set]
        end
    elseif simulation_settings.dc_grid_model == :NonLinear       
        if length(collect(dc_Branch_set)) != 0
            pij_dc = JuMP.value.(model[:pij_dc])
            [grid.DCBranches[branch].PowerFlow_ij_k = [pij_dc[grid.DCBranches[branch].Fr_bus_ID,grid.DCBranches[branch].To_bus_ID,κ] for κ in k] for branch in dc_Branch_set]
            [grid.DCBranches[branch].PowerFlow_ji_k = [pij_dc[grid.DCBranches[branch].To_bus_ID,grid.DCBranches[branch].Fr_bus_ID,κ] for κ in k] for branch in dc_Branch_set]
            [grid.DCBranches[branch].losses_P_k = abs.(grid.DCBranches[branch].PowerFlow_ij_k+grid.DCBranches[branch].PowerFlow_ji_k) for branch in dc_Branch_set]
        end
        if length(prerequisites.dc_Nodes_set) != 0
            v_dc = JuMP.value.(model[:v_dc])
            [grid.DCBuses[bus].V_magnitude_k = [v_dc[bus,κ] for κ in k] for bus in prerequisites.dc_Nodes_set]
        end
    else
        error("Unidentified grid model setting ($simulation_settings.dc_grid_model)")
        return -1
    end

    # Converter data
    if length(prerequisites.converter_set) != 0
        p_conv = JuMP.value.(model[:p_conv])
        [grid.Generators[g].Pg_k = [p_conv[g,κ] for κ in k] for g in collect(union(Set(prerequisites.ac_virtual_gen_set),Set(prerequisites.dc_virtual_gen_set))) ]
    end

    # DC links
    for link_id in keys(grid.DCLinks)
        my_link = grid.DCLinks[link_id]
        my_link_gen_ids = [my_link.Fr_gen_ID,my_link.To_gen_ID]
        p_conv = JuMP.value.(model[:p_conv])
        [grid.Generators[g].Pg_k = [p_conv[g,κ] for κ in k] for g in my_link_gen_ids]
        
        grid.DCLinks[link_id].PowerFlow_ij_k = grid.Generators[my_link.Fr_gen_ID].Pg_k
        grid.DCLinks[link_id].PowerFlow_ji_k = grid.Generators[my_link.To_gen_ID].Pg_k
        grid.DCLinks[link_id].losses_P_k = abs.(grid.Branches[link_id].PowerFlow_ij_k+grid.Branches[link_id].PowerFlow_ji_k)
    end

    # Transmission switching
    if prerequisites.switched_transmission_set != []
        
        z = JuMP.value.(model[:z])
        grid.z_lines = z
        [grid.Branches[branch].GeneralSwitch.SwitchingStatus_k = [z[branch,κ] for κ in k] for branch in prerequisites.switched_transmission_set]
    end

    # Substation reconfiguration and busbar splitting
    if prerequisites.coupler_set != []
        z_l = JuMP.value.(model[:z_l])
        z_c = JuMP.value.(model[:z_c])
        grid.z_reconf = z_l
        grid.z_coupler = z_c
        [grid.Branches[reconf].GeneralSwitch.SwitchingStatus_k = abs.([z_l[reconf,κ] for κ in k]) for reconf in prerequisites.reconf_set]
        [grid.Branches[coupler].GeneralSwitch.SwitchingStatus_k = abs.([z_c[coupler,κ] for κ in k]) for coupler in prerequisites.coupler_set]
        
        for substation_id in keys(grid.Substations)
            if z_c[grid.Substations[substation_id].Reconf_CouplerLines_IDs[1],1] == 1
                grid.Substations[substation_id].is_split = false
            else
                grid.Substations[substation_id].is_split = true
            end
        end

        for bus in prerequisites.aux_bus_set
            if sum(z_l[l,1] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) < 1
                grid.Buses[bus].GeneralSwitch.SwitchingStatus = 0
            end
        end
    end

end

function Generate_Contingencies!(grid::PowerGrid, type::Symbol; include_leafs=false, include_HVDC=false)

    """
    This is a helper function that creates matirces of binary constants to simulate contingencies.
    include_leafs -> whether to include leafs in contingencies or no. If you have leafs in your grid, expect that it will not be N-1 secure
    include_HVDC -> whether to include HVDC lines in contingencies or no
    """

    exception_lines = []

    if !include_leafs
        for bus_id in collect(keys(grid.Buses))
            connected_line_ids = grid.Buses[bus_id].ConnectedLinesIDs
            potential_leaf_id = []
            for branch_id in connected_line_ids
                if grid.Branches[branch_id].BranchType == 0
                    push!(potential_leaf_id, branch_id)
                end
            end
    
            if length(potential_leaf_id) == 1
                push!(exception_lines, potential_leaf_id[1])
            end
    
        end
    end

    N_AC_branches = length([grid.Branches[branch_id] for branch_id in keys(grid.Branches) if grid.Branches[branch_id].BranchType == 0]) - length(exception_lines)
    N_DC_branches = 0
    if include_HVDC
        N_DC_branches = length([grid.DCBranches[branch_id] for branch_id in keys(grid.DCBranches)]) + length([grid.DCLinks[branch_id] for branch_id in keys(grid.DCLinks)])
    end
    N_Line = N_AC_branches + N_DC_branches
    N_Gen = length([grid.Generators[g] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual])
    N_Nodes = grid.N_bus
    N_DC_Nodes = grid.N_dc_bus

    if type == :TC # transmission contingencies usually used with Preventive-SCOPF
        k = N_Line + 1
        a_l = ones(N_Nodes, N_Nodes, k)
        a_l_dc_branch = ones(N_DC_Nodes, N_DC_Nodes, k)
        a_l_dc_link = ones(N_Nodes, N_Nodes, k)
        a_g = ones(N_Gen, k)
        
        c = 1
        [grid.Generators[g].GeneralSwitch.SwitchingStatus_k = ones(k,1) for g in keys(grid.Generators)]
        
        for branch_id in sort(collect(keys(grid.Branches)))
            grid.Branches[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
            if grid.Branches[branch_id].BranchType == 0 && branch_id ∉ exception_lines
                # doing it only for the physical non-auxilliary branches
                c += 1
                grid.Branches[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.Branches[branch_id].Fr_bus_ID
                To_bus_ID = grid.Branches[branch_id].To_bus_ID
                a_l[Fr_bus_ID, To_bus_ID, c] = 0
                a_l[To_bus_ID, Fr_bus_ID, c] = 0
            end
        end

        if include_HVDC
            for branch_id in sort(collect(keys(grid.DCBranches)))
                grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
                c += 1
                grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.DCBranches[branch_id].Fr_bus_ID
                To_bus_ID = grid.DCBranches[branch_id].To_bus_ID
                a_l_dc_branch[Fr_bus_ID, To_bus_ID, c] = 0
                a_l_dc_branch[To_bus_ID, Fr_bus_ID, c] = 0

            end

            for branch_id in sort(collect(keys(grid.DCLinks)))
                grid.DCLinks[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
                c += 1
                grid.DCLinks[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.DCLinks[branch_id].Fr_bus_ID
                To_bus_ID = grid.DCLinks[branch_id].To_bus_ID
                a_l_dc_link[Fr_bus_ID, To_bus_ID, c] = 0
                a_l_dc_link[To_bus_ID, Fr_bus_ID, c] = 0

            end
        end

        if k !== c
            error("Mismatch in contingencies!")
        end

        return a_g, a_l, a_l_dc_branch, a_l_dc_link, k

    elseif type == :TGC # transmission and generation contingencies usually used with Corrective-SCOPF

        k = N_Line + N_Gen + 1
        a_l = ones(N_Nodes, N_Nodes, k)
        a_l_dc_branch = ones(N_DC_Nodes, N_DC_Nodes, k)
        a_l_dc_link = ones(N_Nodes, N_Nodes, k)
        a_g = ones(N_Gen, k)

        c = 1
        
        for branch_id in sort(collect(keys(grid.Branches)))
            grid.Branches[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
            if grid.Branches[branch_id].BranchType == 0 && branch_id ∉ exception_lines
                # doing it only for the physical non-auxilliary branches
                c += 1
                grid.Branches[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.Branches[branch_id].Fr_bus_ID
                To_bus_ID = grid.Branches[branch_id].To_bus_ID
                a_l[Fr_bus_ID, To_bus_ID, c] = 0
                a_l[To_bus_ID, Fr_bus_ID, c] = 0
            end
        end

        if include_HVDC
            for branch_id in sort(collect(keys(grid.DCBranches)))
                grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
                c += 1
                grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.DCBranches[branch_id].Fr_bus_ID
                To_bus_ID = grid.DCBranches[branch_id].To_bus_ID
                a_l_dc_branch[Fr_bus_ID, To_bus_ID, c] = 0
                a_l_dc_branch[To_bus_ID, Fr_bus_ID, c] = 0

            end

            for branch_id in sort(collect(keys(grid.DCLinks)))
                grid.DCLinks[branch_id].GeneralSwitch.SwitchingStatus_k = ones(k,1)
                c += 1
                grid.DCLinks[branch_id].GeneralSwitch.SwitchingStatus_k[c] = 0
                Fr_bus_ID = grid.DCLinks[branch_id].Fr_bus_ID
                To_bus_ID = grid.DCLinks[branch_id].To_bus_ID
                a_l_dc_link[Fr_bus_ID, To_bus_ID, c] = 0
                a_l_dc_link[To_bus_ID, Fr_bus_ID, c] = 0

            end
        end

        for g in keys(grid.Generators)
            grid.Generators[g].GeneralSwitch.SwitchingStatus_k = ones(k,1)
            c += 1
            if grid.Generators[g].GenType  != :virtual
                a_g[g, c] = 0
                grid.Generators[g].GeneralSwitch.SwitchingStatus_k[c] = 0
            end
        end

        if k !== c
            error("Mismatch in contingencies!")
        end

        return a_g, a_l, a_l_dc_branch, a_l_dc_link, k
    end
end

function solve_SCOPF_CCG_model!(grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites; ϵ = 0.01)
    LB = -999999
    UB = 999999

    iter_count = 0
    K_all = prerequisites_data.k
    K_now = [1]
    K_new = []

    t_master = []
    t_sub = []
    LB_it = []
    UB_it = []
    ϵ_it = []
    LS_it = []

    #Main loop
    while abs((UB - LB) / LB) ≥ ϵ
        iter_count = iter_count + 1

        println("Debug msg @ iteration " * string(iter_count))
        # Solve master-problem
        if iter_count == 1
            # prepare system
            System = deepcopy(System_)
            System = split_system!(System, split_buses, "OSR-pre")

            #solve model
            t_master_now = @elapsed LB, z, pg, K_old, LS = solve_master(System, System_, true, [], 1, split_buses)
            println("Length of K_total: ", length(K_old))
            push!(t_master, t_master_now)
            push!(LB_it, LB)
            push!(LS_it, LS)
        else

            # prepare system
            System = deepcopy(System_)
            System = split_system!(System, split_buses, "OSR-OBS")

            println("K_NEW: ", K_new)
            println("K_OLD: ", K_old)

            #solve model
            t_master_now = @elapsed LB, z, pg, K_old, LS = solve_master(System, System_, false, K_old, K_new, split_buses)
            println("Length of K_total: ", length(K_old))
            push!(t_master, t_master_now)
            push!(LB_it, LB)
            push!(LS_it, LS)
        end

        System = deepcopy(System_)
        System = split_system!(System, split_buses, "OSR-OBS")
        # Solve sub-problem
        t_sub_now = @elapsed UB, K_new = solve_sub_problem(System, System_, z, pg, K_old, method, ΔPg, iter_count - 1)
        push!(t_sub, t_sub_now)
        push!(UB_it, UB)
        println("|Iteration Count|LB|UB|")
        println(iter_count, "  ", LB, " ", UB)
    end

end

function SCOPF_MP!(grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites, k_now)

    mod_prereq = deepcopy(prerequisites_data)
    mod_prereq.k = k_now
    model = build_SCOPF_model!(grid, simulation_settings, mod_prereq, optimize=true)

    # LB, z, pg, K_old, LS

end

function SCOPF_SP(grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites, k_prev; n_next_k=1)
    # Threads.@threads for k in prerequisites_data.k
    #     if k != 1
    #         costs[k] = sub_problem(grid, z, pg, a, k)
    #     end
    # end
    # UB = maximum(costs)
    # k_new = _find_next_largest_k(costs, k_prev, n_next_k)
    # return UB, k_new
end

function master_problem!(iteration::Int, grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites; prioritize_splitting=false)
    # Sanity check
    solver = simulation_settings.MILP_solver
    NLP_flag = false
    if simulation_settings.ac_grid_model == :ACOPF || simulation_settings.dc_grid_model == :Nonlinear
        solver = simulation_settings.NLP_solver
        NLP_flag = true
    end

    if NLP_flag
        if simulation_settings.transmission_switching != [] || simulation_settings.substation_switching != []
            error("Switching is not supported with Nonlinear problem settings.")
            return -1
        end
    end

    # Model initialization
    model = Model(solver)
    # model = direct_model(solver())

    # Variables initialization
    single_period_SCOPF_variable_initialization!(model, simulation_settings, prerequisites_data)
    single_period_SCOPF_LoadShedding_variable_initialization!(model, prerequisites_data)

    # AC grid constraints:
    single_period_SC_angle_limits_ac_grid!(prerequisites_data, grid, model)
    single_period_SC_voltage_limits_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_generator_limits_LS_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_nodal_balance_LS_ac_node!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    if prerequisites_data.switched_transmission_set != []
        single_period_SC_switched_transmission_capacity_limits_ac_branch!(prerequisites_data, grid, model)
        single_period_SC_switched_powerflow_ac_branch!(prerequisites_data, grid, model; max_op=simulation_settings.max_transmission_switching)
    end

    single_period_SC_transmission_capacity_limits_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_SC_powerflow_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    if prerequisites_data.Coupler_set != []
        single_period_SC_reconf_split_constraints_ac_grid!(prerequisites_data, grid, model; max_reconf=simulation_settings.max_substation_reconf, max_splitting=simulation_settings.max_busbar_splitting)
    end

    # DC grid constraints:
    if length(keys(grid.DCBuses)) != 0
        single_period_SC_voltage_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_powerflow_dc_branch!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_transmission_capacity_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_SC_nodal_balance_dc_node!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
    end

    # Converter constraints:
    if length(keys(grid.Converters)) != 0 || length(keys(grid.DCLinks)) != 0
        single_period_SC_converter_constraints!(prerequisites_data, grid, model, simulation_settings.converter_model)
    end

    single_period_SC_corrective_measures!(prerequisites_data, grid, model, simulation_settings)

    transmission_switching = prerequisites_data.switched_transmission_set != [] ? true : false
    substation_switching = prerequisites_data.Coupler_set != [] ? true : false
    single_period_SC_LS_objective!(prerequisites_data, grid, model, transmission_switching, substation_switching)

    if prerequisites_data.Coupler_set != [] && prioritize_splitting
        for c in prerequisites_data.Coupler_set
            for k in prerequisites_data.k
                MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), model[:z_c][c, k], 1)
            end
        end
    end
    optimize!(model)

    if !has_values(model)
        error("Problem is INFEASIBLE, iteration -- " * string(iteration))
    end

    SCOPF_post_processing!(grid, model, simulation_settings, prerequisites)
    LB, z, pg, p_conv, K_old, LS = grid.Operating_Cost, grid.z_l, JuMP.value(model[:p][!, 1]), JuMP.value(model[:P_conv][!, 1]), prerequisites_data.k, JuMP.value(model[:p_ls])
    return LB, z, pg, p_conv, K_old, LS
end

function sub_problem(grid::PowerGrid, simulation_settings::SCOPF_SimulationSettings, prerequisites_data::SCOPF_Prerequisites)    
end

function _find_next_largest_k(cost_arr, k_previous, n_next_k)
    continue_itr = true
    found_k = 0
    element_indices = []

    while continue_itr
        a = findfirst(x -> x == maximum(cost_arr), cost_arr)[1]
        if !(a in k_previous)
            push!(element_indices, a)
            arr[a] = 0
            found_k = found_k + 1
        else
            arr[a] = 0
        end
        if found_k == n_next_k || sum(cost_arr) == 0
            continue_itr = false
        end
    end

    return element_indices
end