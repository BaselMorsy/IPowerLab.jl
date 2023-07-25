include("DOPF_Constraints.jl")
include("DOPF_NB.jl")

function Create_Contingency_Maps(grid::PowerGrid, types::Array; include_leafs=false, k=1, converter_modularization=:continuous)

    """
    This is a helper function that creates vectors of binary constants to simulate contingencies.
    `types` -> could be any subset of `[:ac_branch, :dc_branch, :dc_link, :ac_gen, :dc_gen, :conv, :coupler]`
    `include_leafs` -> whether to include leafs in contingencies or no. If you have leafs in your grid, expect that it will not be N-1 secure, but setting this parameter to `false` will disregard these branches
    """

    exception_ac_branches = []
    exception_dc_branches = []

    if !include_leafs
        # this for loop makes sure that no leaf-like AC branch is taken into the contingency considerations
        for bus_id in collect(keys(grid.Buses))
            connected_line_ids = grid.Buses[bus_id].ConnectedLinesIDs
            potential_leaf_id = []
            other_connected_lines = []
            for branch_id in connected_line_ids
                if grid.Branches[branch_id].BranchType == 0
                    push!(potential_leaf_id, branch_id)
                else
                    push!(other_connected_lines, branch_id)
                end
            end

            virtual_gens_at_bus = []
            for gen_id in grid.Buses[bus_id].ConnectedGensIDs
                if grid.Generators[gen_id].GenType == :virtual
                    push!(virtual_gens_at_bus, gen_id)
                end
            end
    
            if length(potential_leaf_id) == 1 && length(virtual_gens_at_bus) == 0 && length(other_connected_lines) == 0 # not connected to any converters or DC links
                push!(exception_ac_branches, potential_leaf_id[1])
            end
    
        end
        exception_ac_branches = collect(Set(exception_ac_branches))
        # this for loop makes sure that no leaf-like DC branch is taken into the contingency considerations
        for dc_bus_id in collect(keys(grid.DCBuses))
            connected_line_ids = grid.DCBuses[dc_bus_id].ConnectedLinesIDs
            potential_leaf_id = []
            for branch_id in connected_line_ids
                push!(potential_leaf_id, branch_id)
            end
            virtual_gens_at_bus = []
            for gen_id in grid.DCBuses[dc_bus_id].ConnectedGensIDs
                if grid.DCGenerators[gen_id].GenType == :virtual
                    push!(virtual_gens_at_bus, gen_id)
                end
            end
    
            if length(potential_leaf_id) == 1 && length(virtual_gens_at_bus) == 0 # not connected to any converters or DC links
                push!(exception_dc_branches, potential_leaf_id[1])
            end
        end
        exception_dc_branches = collect(Set(exception_dc_branches))
    end

    modular_converters = []
    for i in keys(grid.Converter_Duplets)
        for c_id in grid.Converter_Duplets[i]
            push!(modular_converters, c_id)
        end 
    end

    unmodular_converters = []
    for conv_id in sort(collect(keys(grid.Converters)))
        if conv_id ∉ modular_converters
            push!(unmodular_converters, conv_id)
        end
    end

    modular_generators = []
    for i in keys(grid.Generator_Duplets)
        for g_id in grid.Generator_Duplets[i]
            push!(modular_converters, g_id)
        end
    end

    unmodular_generators = []
    for g in keys(grid.Generators)
        if g ∉ modular_generators && grid.Generators[g].GenType != :virtual
            push!(unmodular_generators, g)
        end
    end

    N_AC_branches = length([branch_id for branch_id in keys(grid.Branches) if grid.Branches[branch_id].BranchType == 0]) - length(exception_ac_branches)
    N_DC_branches = length([branch_id for branch_id in keys(grid.DCBranches)]) - length(exception_dc_branches)
    N_DC_Links = length([grid.DCLinks[branch_id] for branch_id in keys(grid.DCLinks)])
    N_AC_gen = length(unmodular_generators) + length(keys(grid.Generator_Duplets))
    N_DC_gen = length([g for g in keys(grid.DCGenerators) if grid.DCGenerators[g].GenType != :virtual])
    if converter_modularization == :discrete
        N_Conv = length([grid.Converters[c] for c in keys(grid.Converters)])
    elseif converter_modularization == :continuous
        N_Conv = length(unmodular_converters) + length(keys(grid.Converter_Duplets))
    end
    N_Couplers = sum([length(grid.Substations[s].Reconf_CouplerLines_IDs) for s in keys(grid.Substations)], init = 0)

    # Sanity check
    @assert length([grid.Converters[c] for c in keys(grid.Converters)]) == length(unmodular_converters) + 2*length(keys(grid.Converter_Duplets))
    @assert length([grid.Generators[g] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]) == length(unmodular_generators) + length(modular_generators)

    N_k = 0
    for contingency_type in types
        if contingency_type == :ac_branch
            N_k += N_AC_branches
        elseif contingency_type == :dc_branch
            N_k += N_DC_branches
        elseif contingency_type == :dc_link
            N_k += N_DC_Links
        elseif contingency_type == :ac_gen
            N_k += N_AC_gen
        elseif contingency_type == :dc_gen
            N_k += N_DC_gen
        elseif contingency_type == :conv
            N_k += N_Conv
        elseif contingency_type == :coupler
            N_k += N_Couplers
        else
            error("No type of contingency matching -> $(contingency_type) !!")
        end
    end
    
    N_e = deepcopy(N_k)

    if k > N_e
        # Sanity check
        error("k ($k) cannot be greater than N ($N_e)")
    end

    if k == 0
        N_k = 1
    else
        N_k = sum(binomial(N_k, k_) for k_ in 1:k) + 1 # covering all contingencies from N-1 to N-k
    end
    
    

    N_AC_branches_ = [branch_id for branch_id in keys(grid.Branches) if grid.Branches[branch_id].BranchType == 0]
    N_DC_branches_ = [branch_id for branch_id in keys(grid.DCBranches)]
    N_DC_Links_    = [branch_id for branch_id in keys(grid.DCLinks)]
    N_AC_gen_      = [g for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]
    N_DC_gen_      = [g for g in keys(grid.DCGenerators) if grid.DCGenerators[g].GenType != :virtual]
    N_conv_        = [c for c in keys(grid.Converters)]
    N_couplers_    = []
    for s in keys(grid.Substations)
        append!(N_couplers_, grid.Substations[s].Reconf_CouplerLines_IDs)
    end
    push!(N_AC_branches_, 0)
    push!(N_DC_branches_, 0)
    push!(N_DC_Links_, 0)
    push!(N_AC_gen_, 0)
    push!(N_DC_gen_, 0)
    push!(N_conv_, 0)
    push!(N_couplers_, 0)

    N_AC_branches_ = maximum(N_AC_branches_)
    N_DC_branches_ = maximum(N_DC_branches_)
    N_DC_Links_    = maximum(N_DC_Links_)
    N_AC_gen_      = maximum(N_AC_gen_)
    N_DC_gen_      = maximum(N_DC_gen_)
    N_conv_        = maximum(N_conv_)
    N_couplers_    = maximum(N_couplers_)

    A_ac_branch = ones(N_AC_branches_, N_k)
    A_dc_branch = ones(N_DC_branches_, N_k)
    A_dc_link   = ones(N_DC_Links_, N_k)
    A_ac_gen    = ones(N_AC_gen_, N_k)
    A_dc_gen    = ones(N_DC_gen_, N_k)
    A_conv      = ones(N_conv_, N_k)
    A_coupler   = ones(N_couplers_, N_k)

    loss_of_generation = []

    if k == 0
        contingency_map = Dict()
        push!(contingency_map, "ac_branch" => A_ac_branch)
        push!(contingency_map, "dc_branch" => A_dc_branch)
        push!(contingency_map, "dc_link" => A_dc_link)
        push!(contingency_map, "ac_gen" => A_ac_gen)
        push!(contingency_map, "dc_gen" => A_dc_gen)
        push!(contingency_map, "conv" => A_conv)
        push!(contingency_map, "coupler" => A_coupler)
        return contingency_map, N_k, loss_of_generation
    end

    c_k = 2
    
    if :ac_branch ∈ types
        # create AC branch contingencies excluding exceptions
        for ac_branch_id in sort(collect(keys(grid.Branches)))
            branch = grid.Branches[ac_branch_id]
            if branch.BranchType == 0 && ac_branch_id ∉ exception_ac_branches
                A_ac_branch[ac_branch_id,c_k] = 0
                c_k += 1
            end
        end
    end

    if :dc_branch ∈ types
        # create DC branch contingencies excluding exceptions
        for dc_branch_id in sort(collect(keys(grid.DCBranches)))
            branch = grid.DCBranches[dc_branch_id]
            if dc_branch_id ∉ exception_dc_branches
                A_dc_branch[dc_branch_id,c_k] = 0
                c_k += 1
            end
        end
    end

    if :dc_link ∈ types
        # create DC link contingencies
        for dc_link_id in sort(collect(keys(grid.DCLinks)))
            dc_link = grid.DCLinks[dc_link_id]
            A_dc_link[dc_link,c_k] = 0
            c_k += 1
        end
    end

    if :ac_gen ∈ types
        # create AC generator contingencies
        for ac_gen_id in sort(unmodular_generators)
            if grid.Generators[ac_gen_id].GenType != :virtual
                A_ac_gen[ac_gen_id,c_k] = 0
                push!(loss_of_generation, c_k)
                c_k += 1
            end
        end
        for g_duplet_id in keys(grid.Generator_Duplets)
            for g_id in grid.Generator_Duplets[g_duplet_id]
                A_conv[g_id,c_k] = 0
            end
            push!(loss_of_generation, c_k)
            c_k += 1 
        end
    end

    if :dc_gen ∈ types
        # create DC generator contingencies
        for dc_gen_id in sort(collect(keys(grid.DCGenerators)))
            if grid.DCGenerators[dc_gen_id].GenType != :virtual
                A_dc_gen[dc_gen_id,c_k] = 0
                push!(loss_of_generation, c_k)
                c_k += 1
            end
        end
    end

    if :conv ∈ types
        # create converter contingencies -> modular converters fall out altogether
        if converter_modularization == :discrete
            for conv_id in sort(collect(keys(grid.Converters)))
                A_conv[conv_id,c_k] = 0
                c_k += 1
            end
        elseif converter_modularization == :continuous
            for conv_id in sort(collect(keys(grid.Converters)))
                if conv_id ∉ modular_converters
                    A_conv[conv_id,c_k] = 0
                    c_k += 1
                end
            end
            for duplet_id in keys(grid.Converter_Duplets)
                for c_id in grid.Converter_Duplets[duplet_id]
                    A_conv[c_id,c_k] = 0
                end
                c_k += 1 
            end
        end
    end

    if :coupler ∈ types
        # create coupler contingencies
        for substation_id in sort(collect(keys(grid.Substations)))
            coupler_ids = grid.Substations[substation_id].Reconf_CouplerLines_IDs
            for coupler_id in coupler_ids
                A_coupler[coupler_id,c_k] = 0
                c_k += 1
            end
        end
    end

    if k == 1

        @assert c_k == N_k + 1
        
        contingency_map = Dict()
        push!(contingency_map, "ac_branch" => A_ac_branch)
        push!(contingency_map, "dc_branch" => A_dc_branch)
        push!(contingency_map, "dc_link" => A_dc_link)
        push!(contingency_map, "ac_gen" => A_ac_gen)
        push!(contingency_map, "dc_gen" => A_dc_gen)
        push!(contingency_map, "conv" => A_conv)
        push!(contingency_map, "coupler" => A_coupler)

        return contingency_map, N_k, loss_of_generation
    else

        Contingent_Elements = collect(2:N_e)

        for k_ in collect(2:k)
            contingency_sets = collect(combinations(Contingent_Elements,k_))
            for contingency_set in contingency_sets

                for ind in contingency_set
                    if ind in loss_of_generation
                        push!(loss_of_generation, ind)
                        break
                    end
                end

                vectors = [A_ac_branch[:,i] for i in contingency_set]
                A_ac_branch[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_dc_branch[:,i] for i in contingency_set]
                A_dc_branch[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_dc_link[:,i] for i in contingency_set]
                A_dc_link[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_ac_gen[:,i] for i in contingency_set]
                A_ac_gen[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_dc_gen[:,i] for i in contingency_set]
                A_dc_gen[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_conv[:,i] for i in contingency_set]
                A_conv[:,c_k] = [minimum([values...]) for values in zip(vectors...)]


                vectors = [A_coupler[:,i] for i in contingency_set]
                A_coupler[:,c_k] = [minimum([values...]) for values in zip(vectors...)]
                
                c_k += 1 # move to new contingency set
            end
        end

        @assert c_k == N_k + 1
        
        contingency_map = Dict()
        push!(contingency_map, "ac_branch" => A_ac_branch)
        push!(contingency_map, "dc_branch" => A_dc_branch)
        push!(contingency_map, "dc_link" => A_dc_link)
        push!(contingency_map, "ac_gen" => A_ac_gen)
        push!(contingency_map, "dc_gen" => A_dc_gen)
        push!(contingency_map, "conv" => A_conv)
        push!(contingency_map, "coupler" => A_coupler)

        return contingency_map, N_k, loss_of_generation
    end
end

function compile_grid_topology!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, switching_specs::Dict;  modular_converters=[])
    for conv_id in modular_converters
        split_element!(grid, :Converter, conv_id; n=2, modularization=SimulationSettings.converter_modularization)
    end

    # Preparing the grid topology by adding needed auxilliary nodes and branches
    substation_ids = switching_specs["substation_ids"]
    hybrid_substation_ids = switching_specs["hybrid_substation_ids"]
    B2B_cap = switching_specs["B2B_cap"]

    convert_bus2substation!(grid,substation_ids,false) # Regular substations
    convert_bus2substation!(grid,hybrid_substation_ids,true,B2B_cap) # Soft switch substations
end

function preprocess_bids!(grid::PowerGrid, original_order_book::OrderBook, SimulationSettings::DOPF_SimulationSettings)
    order_book = deepcopy(original_order_book)
    ac_gen_id_to_gen_root = Dict()
    root_gen_to_duplicate_gen = Dict() # useful for deciding on commitable and uncommitable generators because the duplicate_gen follows the root_gen
    all_gens = deepcopy(sort(collect(keys(order_book.Gen_bids))))
    for gen_id in all_gens
        push!(ac_gen_id_to_gen_root, gen_id => gen_id)
        some_list = [order_book.Gen_bids[gen_id]["price"][t] for t in SimulationSettings.time_horizon]
        if any(x -> x > 1, map(length, some_list) )
            # assuming that two prices only happen in case of down-regulation (a generator will be required to decrease their output)
            split_element!(grid, :Generator, gen_id)
            new_gen_id = grid.N_gen
            push!(ac_gen_id_to_gen_root, new_gen_id => gen_id)
            push!(root_gen_to_duplicate_gen, gen_id => new_gen_id)
            new_gen_bid_price = Dict()
            new_gen_bid_qty = Dict()
            for t in SimulationSettings.time_horizon
                down_reg_price = deepcopy(order_book.Gen_bids[gen_id]["price"][t][2])
                down_reg_qty = deepcopy(order_book.Gen_bids[gen_id]["qty"][t][2])
                # up-reg generator
                order_book.Gen_bids[gen_id]["price"][t] = [deepcopy(order_book.Gen_bids[gen_id]["price"][t][1])]
                order_book.Gen_bids[gen_id]["qty"][t][2] = 0

                # down-reg generator
                push!(new_gen_bid_price, t => [down_reg_price])
                push!(new_gen_bid_qty, t => [0, down_reg_qty])
            end
            new_gen_bid = Dict("price" => new_gen_bid_price, "qty" => new_gen_bid_qty)
            push!(order_book.Gen_bids, new_gen_id => new_gen_bid)
            
        end
    end

    order_book.Gen_ids = [g for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]
    initialize_session!(order_book)
    return ac_gen_id_to_gen_root, root_gen_to_duplicate_gen, order_book
end

function compile_prerequisites_DOPF!(grid::PowerGrid, market::PowerMarket, order_book::OrderBook, SimulationSettings::DOPF_SimulationSettings, contingency_specs::Dict, switching_specs::Dict,
        generation_specs::Dict; relaxed_physics_lines=[], relaxed_physics_nodes=[], relaxed_capacity_lines=[],reference_node=nothing)
    
    contingency_map, N_k, loss_of_generation = Create_Contingency_Maps(grid, SimulationSettings.contingency_types; include_leafs=contingency_specs["include_leafs"], k=contingency_specs["k"])
    K = sort(collect(1:N_k))


    ac_transmission_branches = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 0]
    ac_all_reconf_line_ids = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 1]
    ac_all_coupler_line_ids = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 2]
    
    ac_node_ids = collect(keys(grid.Buses))
    ac_aux_bus_ids = [key for key in ac_node_ids if grid.Buses[key].BusType == 1]
    ac_branch_ids = collect(keys(grid.Branches))
    ac_gen_ids = [g for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]
    ac_load_ids = collect(keys(grid.Loads))
    ac_load_shedding_ids = ac_load_ids
    ac_gen_id_to_gen_root = generation_specs["ac_gen_id_to_gen_root"]
    root_gen_to_duplicate_gen = generation_specs["root_gen_to_duplicate_gen"]

    if SimulationSettings.contingency_redispatch_condition == :all
        contingency_redispatch = K
    elseif SimulationSettings.contingency_redispatch_condition == :loss_of_generation
        contingency_redispatch = loss_of_generation
    elseif SimulationSettings.contingency_redispatch_condition == :none
        contingency_redispatch = []
    else
        error("Invalid condition for contingency redispatch: $(SimulationSettings.contingency_redispatch_condition)")
    end

    commitable_gen_ids = generation_specs["commitable_gen_ids"]
    non_commitable_gen_ids = generation_specs["non_commitable_gen_ids"]
    fixed_commitments = generation_specs["fixed_commitments"]
    fixed_schedules = generation_specs["fixed_schedules"]
    
    #### Sanity check for commitment variables ####
    empty_set = []
    append!(empty_set, commitable_gen_ids)
    append!(empty_set, non_commitable_gen_ids)
    append!(empty_set, collect(keys(fixed_commitments)))
    append!(empty_set, collect(keys(fixed_schedules)))
    @assert union(Set(commitable_gen_ids),Set(non_commitable_gen_ids),Set(collect(keys(fixed_commitments))),Set(collect(keys(fixed_schedules)))) == Set(ac_gen_ids)
    @assert length(empty_set) == length(commitable_gen_ids) + length(non_commitable_gen_ids) + length(collect(keys(fixed_commitments))) + length(collect(keys(fixed_schedules)))
    ################################################

    ac_switchable_branch_set = union(Set(switching_specs["ac_active_dynamic_branch_ids"]), Set(switching_specs["ac_fixed_dynamic_branch_ids"]))
    ac_static_branch_ids = collect(setdiff(Set(ac_transmission_branches),ac_switchable_branch_set))
    ac_active_dynamic_branch_ids = switching_specs["ac_active_dynamic_branch_ids"]
    ac_active_reconf_ids = ac_all_reconf_line_ids
    ac_active_coupler_ids = ac_all_coupler_line_ids
    ac_fixed_dynamic_branch_ids = switching_specs["ac_fixed_dynamic_branch_ids"]
    ac_fixed_reconf_ids = []
    ac_fixed_coupler_ids = []
    fixed_topology = switching_specs["fixed_topology"]
    normally_opened_reconf_lines = []
    for substation_id in keys(grid.Substations)
        section_ids = grid.Substations[substation_id].BusbarSections_IDs
        for section_id in section_ids[2:end]
            connected_lines = grid.Buses[section_id].ConnectedLinesIDs
            for line_id in connected_lines
                push!(normally_opened_reconf_lines,line_id)
            end
        end
    end

    dc_link_ids = collect(keys(grid.DCLinks))
    dc_link_vritual_gen_ids = []
    for dc_link_id in dc_link_ids
        push!(dc_link_vritual_gen_ids, grid.DCLinks[dc_link_id].Fr_gen_ID)
        push!(dc_link_vritual_gen_ids, grid.DCLinks[dc_link_id].To_gen_ID)
    end
    conv_ids = [conv_id for conv_id in keys(grid.Converters) if grid.Converters[conv_id].type == :ACDC]
    b2b_conv_ids = [conv_id for conv_id in keys(grid.Converters) if grid.Converters[conv_id].type == :B2B]
    b2b_gen_ids = []
    for b2b_id in b2b_conv_ids
        push!(b2b_gen_ids, grid.Converters[b2b_id].gen_dc_id)
        push!(b2b_gen_ids, grid.Converters[b2b_id].gen_ac_id)
    end
    conv_ac_side_virtual_gen_ids = []
    conv_dc_side_virtual_gen_ids = []
    ac_gen_to_conv_id_dict = Dict()
    dc_gen_to_conv_id_dict = Dict()
    for conv_id in conv_ids
        push!(conv_ac_side_virtual_gen_ids, grid.Converters[conv_id].gen_ac_id)
        push!(conv_dc_side_virtual_gen_ids, grid.Converters[conv_id].gen_dc_id)
        push!(dc_gen_to_conv_id_dict, grid.Converters[conv_id].gen_dc_id => conv_id)
        push!(ac_gen_to_conv_id_dict, grid.Converters[conv_id].gen_ac_id => conv_id)
    end
    conv_duplets = grid.Converter_Duplets # dictionary of converter duplets


    dc_node_ids = collect(keys(grid.DCBuses))
    dc_branch_ids = collect(keys(grid.DCBranches))
    dc_gen_ids = [g for g in keys(grid.DCGenerators) if grid.DCGenerators[g].GenType != :virtual]
    dc_load_ids = collect(keys(grid.DCLoads))

    Schedule = market.Schedule
    Order_Book = order_book

    Contingency_Map = contingency_map
    k = K # vector of contingency indices
    k_t = Dict() # dictionary of contingencies considered at time t
    base_contingency = [1]
    for t in SimulationSettings.time_horizon
        if SimulationSettings.Meta_solver == :CCG
            push!(k_t, t => append!(deepcopy(base_contingency), loss_of_generation))
        elseif SimulationSettings.Meta_solver == :none
            push!(k_t, t => k)
        end
    end

    M_l = Dict()
    for l in ac_switchable_branch_set
        push!(M_l, l => 1.2*grid.S_base*(1/grid.Branches[l].x))
    end
    M_δ = 1.2
    M_E = grid.S_base*10

    ####

    return DOPF_Prerequisites(time_horizon=SimulationSettings.time_horizon, base_MVA=grid.S_base, ac_node_ids=ac_node_ids, ac_aux_bus_ids=ac_aux_bus_ids, ac_branch_ids=ac_branch_ids,
        ac_gen_ids=ac_gen_ids, ac_load_ids=ac_load_ids, ac_load_shedding_ids=ac_load_shedding_ids,
        ac_gen_id_to_gen_root=ac_gen_id_to_gen_root, root_gen_to_duplicate_gen=root_gen_to_duplicate_gen, contingency_redispatch=contingency_redispatch,
        commitable_gen_ids=commitable_gen_ids, non_commitable_gen_ids=non_commitable_gen_ids, fixed_commitments=fixed_commitments,
        fixed_schedules=fixed_schedules, ac_static_branch_ids=ac_static_branch_ids, ac_active_dynamic_branch_ids=ac_active_dynamic_branch_ids, ac_active_reconf_ids=ac_active_reconf_ids,
        ac_active_coupler_ids=ac_active_coupler_ids, ac_fixed_dynamic_branch_ids=ac_fixed_dynamic_branch_ids, ac_fixed_reconf_ids=ac_fixed_reconf_ids, ac_fixed_coupler_ids=ac_fixed_coupler_ids,
        fixed_topology=fixed_topology, normally_opened_reconf_lines=normally_opened_reconf_lines,
        dc_link_ids=dc_link_ids, dc_link_vritual_gen_ids=dc_link_vritual_gen_ids, conv_ids=conv_ids, b2b_conv_ids=b2b_conv_ids, b2b_gen_ids=b2b_gen_ids,
        conv_ac_side_virtual_gen_ids=conv_ac_side_virtual_gen_ids,conv_dc_side_virtual_gen_ids=conv_dc_side_virtual_gen_ids, ac_gen_to_conv_id_dict=ac_gen_to_conv_id_dict, dc_gen_to_conv_id_dict=dc_gen_to_conv_id_dict,
        conv_duplets=conv_duplets, dc_node_ids=dc_node_ids, dc_branch_ids=dc_branch_ids, dc_gen_ids=dc_gen_ids, dc_load_ids=dc_load_ids, Schedule=Schedule, Order_Book=Order_Book,
        Contingency_Map=Contingency_Map, k=k, k_t=k_t, M_l=M_l, M_δ=M_δ, M_E=M_E, relaxed_physics_lines=relaxed_physics_lines, relaxed_physics_nodes=relaxed_physics_nodes,
        relaxed_capacity_lines=relaxed_capacity_lines, reference_node=reference_node)

end

function compile_prerequisites_DOPF_no_market!(grid::PowerGrid, simulation_type::Symbol, SimulationSettings::DOPF_SimulationSettings, contingency_specs::Dict, switching_specs::Dict,
    generation_specs::Dict; relaxed_physics_lines=[], relaxed_physics_nodes=[], relaxed_capacity_lines=[],reference_node=nothing)

    """
    simulation_type: could be => :ED - :OPF - :SCOPF - :UC - :NCUC - :SCUC
    """
    gen_ids = [g for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]
    load_ids = collect(keys(grid.Loads))

    # create a market with base grid
    market = create_market(gen_ids, load_ids; name="MyMarket", horizon=SimulationSettings.time_horizon)
    attach_grid2market!(market, deepcopy(grid))

    # create an order book and add bids automatically from grid
    order_book = OrderBook(Market="day-ahead", Clearing=simulation_type, Pricing=:Uniform, horizon=SimulationSettings.time_horizon, Gen_ids=gen_ids, Load_ids=load_ids)
    initialize_session!(order_book)
    create_auto_base_bids!(market, order_book; new_day=true)

    ac_gen_id_to_gen_root, root_gen_to_duplicate_gen, order_book_new = preprocess_bids!(grid, order_book, SimulationSettings)
    all_gen_ids = collect(keys(ac_gen_id_to_gen_root))
    
    contingency_map, N_k, loss_of_generation = Create_Contingency_Maps(grid, SimulationSettings.contingency_types; include_leafs=contingency_specs["include_leafs"], k=contingency_specs["k"])
    K = sort(collect(1:N_k))


    ac_transmission_branches = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 0]
    ac_all_reconf_line_ids = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 1]
    ac_all_coupler_line_ids = [key for key in keys(grid.Branches) if grid.Branches[key].BranchType == 2]
    
    ac_node_ids = collect(keys(grid.Buses))
    ac_aux_bus_ids = [key for key in ac_node_ids if grid.Buses[key].BusType == 1]
    ac_branch_ids = collect(keys(grid.Branches))
    ac_gen_ids = [g for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual]
    ac_load_ids = collect(keys(grid.Loads))
    ac_load_shedding_ids = ac_load_ids

    if SimulationSettings.contingency_redispatch_condition == :all
        contingency_redispatch = K
    elseif SimulationSettings.contingency_redispatch_condition == :loss_of_generation
        contingency_redispatch = loss_of_generation
    elseif SimulationSettings.contingency_redispatch_condition == :none
        contingency_redispatch = []
    else
        error("Invalid condition for contingency redispatch: $(SimulationSettings.contingency_redispatch_condition)")
    end

    commitable_gen_ids = generation_specs["commitable_gen_ids"]
    non_commitable_gen_ids = generation_specs["non_commitable_gen_ids"]
    fixed_commitments = generation_specs["fixed_commitments"]
    fixed_schedules = generation_specs["fixed_schedules"]
    
    #### Sanity check for commitment variables ####
    empty_set = []
    append!(empty_set, commitable_gen_ids)
    append!(empty_set, non_commitable_gen_ids)
    append!(empty_set, collect(keys(fixed_commitments)))
    append!(empty_set, collect(keys(fixed_schedules)))
    @assert union(Set(commitable_gen_ids),Set(non_commitable_gen_ids),Set(collect(keys(fixed_commitments))),Set(collect(keys(fixed_schedules)))) == Set(ac_gen_ids)
    @assert length(empty_set) == length(commitable_gen_ids) + length(non_commitable_gen_ids) + length(collect(keys(fixed_commitments))) + length(collect(keys(fixed_schedules)))
    ################################################

    ac_switchable_branch_set = union(Set(switching_specs["ac_active_dynamic_branch_ids"]), Set(switching_specs["ac_fixed_dynamic_branch_ids"]))
    ac_static_branch_ids = collect(setdiff(Set(ac_transmission_branches),ac_switchable_branch_set))
    ac_active_dynamic_branch_ids = switching_specs["ac_active_dynamic_branch_ids"]
    ac_active_reconf_ids = ac_all_reconf_line_ids
    ac_active_coupler_ids = ac_all_coupler_line_ids
    ac_fixed_dynamic_branch_ids = switching_specs["ac_fixed_dynamic_branch_ids"]
    ac_fixed_reconf_ids = []
    ac_fixed_coupler_ids = []
    fixed_topology = switching_specs["fixed_topology"]
    normally_opened_reconf_lines = []
    for substation_id in keys(grid.Substations)
        section_ids = grid.Substations[substation_id].BusbarSections_IDs
        for section_id in section_ids[2:end]
            connected_lines = grid.Buses[section_id].ConnectedLinesIDs
            for line_id in connected_lines
                push!(normally_opened_reconf_lines,line_id)
            end
        end
    end

    dc_link_ids = collect(keys(grid.DCLinks))
    dc_link_vritual_gen_ids = []
    for dc_link_id in dc_link_ids
        push!(dc_link_vritual_gen_ids, grid.DCLinks[dc_link_id].Fr_gen_ID)
        push!(dc_link_vritual_gen_ids, grid.DCLinks[dc_link_id].To_gen_ID)
    end
    conv_ids = [conv_id for conv_id in keys(grid.Converters) if grid.Converters[conv_id].type == :ACDC]
    b2b_conv_ids = [conv_id for conv_id in keys(grid.Converters) if grid.Converters[conv_id].type == :B2B]
    b2b_gen_ids = []
    for b2b_id in b2b_conv_ids
        push!(b2b_gen_ids, grid.Converters[b2b_id].gen_dc_id)
        push!(b2b_gen_ids, grid.Converters[b2b_id].gen_ac_id)
    end
    conv_ac_side_virtual_gen_ids = []
    conv_dc_side_virtual_gen_ids = []
    ac_gen_to_conv_id_dict = Dict()
    dc_gen_to_conv_id_dict = Dict()
    for conv_id in conv_ids
        push!(conv_ac_side_virtual_gen_ids, grid.Converters[conv_id].gen_ac_id)
        push!(conv_dc_side_virtual_gen_ids, grid.Converters[conv_id].gen_dc_id)
        push!(dc_gen_to_conv_id_dict, grid.Converters[conv_id].gen_dc_id => conv_id)
        push!(ac_gen_to_conv_id_dict, grid.Converters[conv_id].gen_ac_id => conv_id)
    end
    conv_duplets = grid.Converter_Duplets # dictionary of converter duplets


    dc_node_ids = collect(keys(grid.DCBuses))
    dc_branch_ids = collect(keys(grid.DCBranches))
    dc_gen_ids = [g for g in keys(grid.DCGenerators) if grid.DCGenerators[g].GenType != :virtual]
    dc_load_ids = collect(keys(grid.DCLoads))

    Schedule = market.Schedule
    Order_Book = order_book

    Contingency_Map = contingency_map
    k = K # vector of contingency indices
    k_t = Dict() # dictionary of contingencies considered at time t
    base_contingency = [1]
    for t in SimulationSettings.time_horizon
        if SimulationSettings.Meta_solver == :CCG
            push!(k_t, t => append!(deepcopy(base_contingency), loss_of_generation))
        elseif SimulationSettings.Meta_solver == :none
            push!(k_t, t => k)
        end
    end

    M_l = Dict()
    for l in ac_switchable_branch_set
        push!(M_l, l => 1.2*grid.S_base*(1/grid.Branches[l].x))
    end
    M_δ = 1.2
    M_E = grid.S_base*10

    ####

    return DOPF_Prerequisites(time_horizon=SimulationSettings.time_horizon, base_MVA=grid.S_base, ac_node_ids=ac_node_ids, ac_aux_bus_ids=ac_aux_bus_ids, ac_branch_ids=ac_branch_ids,
        ac_gen_ids=ac_gen_ids, ac_load_ids=ac_load_ids, ac_load_shedding_ids=ac_load_shedding_ids,
        ac_gen_id_to_gen_root=ac_gen_id_to_gen_root, root_gen_to_duplicate_gen=root_gen_to_duplicate_gen, contingency_redispatch=contingency_redispatch,
        commitable_gen_ids=commitable_gen_ids, non_commitable_gen_ids=non_commitable_gen_ids, fixed_commitments=fixed_commitments,
        fixed_schedules=fixed_schedules, ac_static_branch_ids=ac_static_branch_ids, ac_active_dynamic_branch_ids=ac_active_dynamic_branch_ids, ac_active_reconf_ids=ac_active_reconf_ids,
        ac_active_coupler_ids=ac_active_coupler_ids, ac_fixed_dynamic_branch_ids=ac_fixed_dynamic_branch_ids, ac_fixed_reconf_ids=ac_fixed_reconf_ids, ac_fixed_coupler_ids=ac_fixed_coupler_ids,
        fixed_topology=fixed_topology, normally_opened_reconf_lines=normally_opened_reconf_lines,
        dc_link_ids=dc_link_ids, dc_link_vritual_gen_ids=dc_link_vritual_gen_ids, conv_ids=conv_ids, b2b_conv_ids=b2b_conv_ids, b2b_gen_ids=b2b_gen_ids,
        conv_ac_side_virtual_gen_ids=conv_ac_side_virtual_gen_ids,conv_dc_side_virtual_gen_ids=conv_dc_side_virtual_gen_ids, ac_gen_to_conv_id_dict=ac_gen_to_conv_id_dict, dc_gen_to_conv_id_dict=dc_gen_to_conv_id_dict,
        conv_duplets=conv_duplets, dc_node_ids=dc_node_ids, dc_branch_ids=dc_branch_ids, dc_gen_ids=dc_gen_ids, dc_load_ids=dc_load_ids, Schedule=Schedule, Order_Book=Order_Book,
        Contingency_Map=Contingency_Map, k=k, k_t=k_t, M_l=M_l, M_δ=M_δ, M_E=M_E, relaxed_physics_lines=relaxed_physics_lines, relaxed_physics_nodes=relaxed_physics_nodes,
        relaxed_capacity_lines=relaxed_capacity_lines, reference_node=reference_node), market
end

function run_DOPF_simulation!(grid::PowerGrid, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites, virtual_market::PowerMarket; update_grid=true)
    return clear_market!(virtual_market, grid, deepcopy(prerequisites_data.Order_Book), simulation_settings,
        prerequisites_data; update_grid=update_grid)
end

function DOPF_post_processing!(model::Model, grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites, order_book::OrderBook)
    
    Branch_set = keys(grid.Branches)
    dc_Branch_set = keys(grid.DCBranches)
    k = sort(collect(prerequisites_data.k))

    if haskey(model, :δ)
        for bus_id in prerequisites_data.ac_node_ids
            if grid.Buses[bus_id].δ_tk == []
                grid.Buses[bus_id].δ_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Buses[bus_id].δ_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Buses[bus_id].δ_tk[t], k => JuMP.value.(model[:δ][bus_id,k,t]))
                end
            end
        end    
    end

    if haskey(model, :p_branch_ac)
        for branch_id in prerequisites_data.ac_branch_ids
            if grid.Branches[branch_id].PowerFlow_ij_tk == []
                grid.Branches[branch_id].PowerFlow_ij_tk = Dict()
                grid.Branches[branch_id].PowerFlow_ji_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[branch_id].PowerFlow_ij_tk, t => Dict())
                push!(grid.Branches[branch_id].PowerFlow_ji_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[branch_id].PowerFlow_ij_tk[t], k => JuMP.value.(model[:p_branch_ac][branch_id,1,k,t]))
                    push!(grid.Branches[branch_id].PowerFlow_ji_tk[t], k => JuMP.value.(model[:p_branch_ac][branch_id,2,k,t]))
                end
            end
        end  
    end

    if haskey(model, :p_gen_ac)
        for ac_gen_id in prerequisites_data.ac_gen_ids
            
            grid.Generators[ac_gen_id].Pg_tk = Dict()
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[ac_gen_id].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Generators[ac_gen_id].Pg_tk[t], k => JuMP.value.(model[:p_gen_ac][ac_gen_id,k,t]) + get(get(prerequisites_data.Schedule["gen"],ac_gen_id,Dict()),t,0))
                end
            end

        end
    end

    if haskey(model, :p_ls_ac)
        for ac_ls_id in prerequisites_data.ac_load_shedding_ids
            if grid.Loads[ac_ls_id].P_shedding_tk == []
                grid.Loads[ac_ls_id].P_shedding_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Loads[ac_ls_id].P_shedding_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Loads[ac_ls_id].P_shedding_tk[t], k => JuMP.value.(model[:p_ls_ac][ac_ls_id,k,t]))
                end
            end
        end
    end

    if haskey(model, :u_gt)
        for g in prerequisites_data.commitable_gen_ids
            if grid.Generators[g].u_gt == []
                grid.Generators[g].u_gt = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[g].u_gt, t => JuMP.value.(model[:u_gt][g, t]))
            end
        end
    end

    if haskey(model, :u_gt_f)
        for g in keys(prerequisites_data.fixed_commitments)
            if grid.Generators[g].u_gt == []
                grid.Generators[g].u_gt = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[g].u_gt, t => JuMP.value.(model[:u_gt_f][g, t]))
            end
        end
    end

    if haskey(model, :z_l)
        for l in prerequisites_data.ac_active_dynamic_branch_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_l][l,k,t]))
                end
            end
        end
        grid.z_lines = JuMP.value.(model[:z_l])
    end

    if haskey(model, :z_l_f)
        for l in prerequisites_data.ac_fixed_dynamic_branch_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_l_f][l,k,t]))
                end
            end
        end
    end

    if haskey(model, :z_r)
        for l in prerequisites_data.ac_active_reconf_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_r][l,k,t]))
                end
            end
        end
        grid.z_reconf = JuMP.value.(model[:z_r])
    end

    if haskey(model, :z_c)
        for l in prerequisites_data.ac_active_coupler_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_c][l,k,t]))
                end
            end
        end
        grid.z_coupler = JuMP.value.(model[:z_c])
    end

    if haskey(model, :z_r_f)
        for l in prerequisites_data.ac_fixed_reconf_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_r_f][l,k,t]))
                end
            end
        end
        
    end

    if haskey(model, :z_c_f)
        for l in prerequisites_data.ac_fixed_coupler_ids
            if grid.Branches[l].GeneralSwitch.SwitchingStatus_tk == []
                grid.Branches[l].GeneralSwitch.SwitchingStatus_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Branches[l].GeneralSwitch.SwitchingStatus_tk[t], k => JuMP.value.(model[:z_c_f][l,k,t]))
                end
            end
        end
    end
 
    if haskey(model, :p_branch_dc)
        for branch_id in prerequisites_data.dc_branch_ids
            if grid.DCBranches[branch_id].PowerFlow_ij_tk == []
                grid.DCBranches[branch_id].PowerFlow_ij_tk = Dict()
                grid.DCBranches[branch_id].PowerFlow_ji_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.DCBranches[branch_id].PowerFlow_ij_tk, t => Dict())
                push!(grid.DCBranches[branch_id].PowerFlow_ji_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.DCBranches[branch_id].PowerFlow_ij_tk[t], k => JuMP.value.(model[:p_branch_dc][branch_id,1,k,t]))
                    push!(grid.DCBranches[branch_id].PowerFlow_ji_tk[t], k => JuMP.value.(model[:p_branch_dc][branch_id,2,k,t]))
                end
            end
        end  
    end

    if haskey(model, :p_gen_dc)
        for dc_gen_id in prerequisites_data.dc_gen_ids
            if grid.DCGenerators[g].Pg_tk == []
                grid.DCGenerators[dc_gen_id].Pg_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.DCGenerators[dc_gen_id].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.DCGenerators[dc_gen_id].Pg_tk[t], k => JuMP.value.(model[:p_gen_ac][dc_gen_id,k,t]))
                end
            end
        end
    end

    if haskey(model, :p_conv_dc)
        for g in prerequisites_data.conv_dc_side_virtual_gen_ids
            if grid.DCGenerators[g].Pg_tk == []
                grid.DCGenerators[g].Pg_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.DCGenerators[g].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.DCGenerators[g].Pg_tk[t], k => JuMP.value.(model[:p_conv_dc][g,k,t]))
                end
            end
        end
    end

    if haskey(model, :p_conv_ac)
        for g in prerequisites_data.conv_ac_side_virtual_gen_ids
            if grid.Generators[g].Pg_tk == []
                grid.Generators[g].Pg_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[g].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Generators[g].Pg_tk[t], k => JuMP.value.(model[:p_conv_ac][g,k,t]))
                end
            end
        end
    end

    if haskey(model, :p_conv_b2b)
        for g in prerequisites_data.b2b_gen_ids
            if grid.Generators[g].Pg_tk == []
                grid.Generators[g].Pg_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[g].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Generators[g].Pg_tk[t], k => JuMP.value.(model[:p_conv_b2b][g,k,t]))
                end
            end
        end
    end

    if haskey(model, :γ_conv)
        if ! haskey(grid.misc, "converter_allocation_by_duplete_id")
            push!(grid.misc, "converter_allocation_by_duplete_id" => Dict())
        end

        for duplet_id in collect(keys(prerequisites_data.conv_duplets))
            if ! haskey(grid.misc["converter_allocation_by_duplete_id"], duplet_id)
                push!(grid.misc["converter_allocation_by_duplete_id"], duplet_id => Dict())
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.misc["converter_allocation_by_duplete_id"][duplet_id], t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.misc["converter_allocation_by_duplete_id"][duplet_id][t], k => JuMP.value.(model[:γ_conv][duplet_id,k,t]))
                end
            end
        end
    end

    if haskey(model, :p_dclink)
        for g in prerequisites_data.dc_link_vritual_gen_ids
            if grid.Generators[g].Pg_tk == []
                grid.Generators[g].Pg_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.Generators[g].Pg_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.Generators[g].Pg_tk[t], k => JuMP.value.(model[:p_dclink][g,k,t]))
                end
            end
        end

        for dc_link in prerequisites_data.dc_link_ids
            g_fr = grid.DCLinks[dc_link].Fr_gen_ID
            g_to = grid.DCLinks[dc_link].To_gen_ID
            if grid.DCLinks[dc_link].PowerFlow_ij_tk == []
                grid.DCLinks[dc_link].PowerFlow_ij_tk = Dict()
                grid.DCLinks[dc_link].PowerFlow_ji_tk = Dict()
            end
            for t in keys(prerequisites_data.k_t)
                push!(grid.DCLinks[dc_link].PowerFlow_ij_tk, t => Dict())
                push!(grid.DCLinks[dc_link].PowerFlow_ji_tk, t => Dict())
                for k in prerequisites_data.k_t[t]
                    push!(grid.DCLinks[dc_link].PowerFlow_ij_tk[t], k => grid.Generators[g_fr].Pg_tk[t][k])
                    push!(grid.DCLinks[dc_link].PowerFlow_ji_tk[t], k => grid.Generators[g_to].Pg_tk[t][k])
                end
            end
        end
    end

end

function update_grid_tables_DOPF!(grid::PowerGrid; t=1, k=1)
    grid.BusData_output = DataFrame(Bus = Int64[], V = Float64[], δ = Float64[],
        Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[], P_ls = Float64[], Q_ls = Float64[],type=[])
    
    grid.LineLoading = DataFrame(BranchID = Int64[], FromBus = Int64[],ToBus = Int64[]
        ,PL_fr_to = Float64[],PL_to_fr = Float64[],
        PLoss = Float64[],QL_fr_to = Float64[]
        ,QL_to_fr = Float64[],QLoss = Float64[], Utilization = Float64[],type=[])

    grid.Converter_flow = DataFrame(ConverterID = Int64[], AC_Bus = Int64[], DC_Bus = Int64[],
        P_ACDC = Float64[], P_DCAC = Float64[], PLoss = Float64[], Utilization = Float64[],type = [])

    # Populate "BusData_output" dataframe
    # Bus = Int64[], V = Float64[], δ = Float64[],Pg = Float64[],Qg = Float64[],Pd = Float64[], Qd = Float64[],P_ls = Float64[], Q_ls = Float64[]
    for bus_id in sort(collect(keys(grid.Buses)))
        my_bus = grid.Buses[bus_id]

        V = get(get(my_bus.V_magnitude_tk,t,Dict()),k,1)
        δ = get(get(my_bus.δ_tk,t,Dict()),k,-1)
        Pg = sum([get(get(grid.Generators[g].Pg_tk,t,Dict()),k,0) for g in my_bus.ConnectedGensIDs if grid.Generators[g].GenType != :virtual], init=0)
        Qg = sum([get(get(grid.Generators[g].Qg_tk,t,Dict()),k,0) for g in my_bus.ConnectedGensIDs if grid.Generators[g].GenType != :virtual], init=0)
        Pd = sum([grid.Loads[d].Pd_t[t] for d in my_bus.ConnectedLoadsIDs], init=0)
        Qd = sum([grid.Loads[d].Qd_t[t] for d in my_bus.ConnectedLoadsIDs], init=0)
        P_ls = sum([get(get(grid.Loads[d].P_shedding_tk,t,Dict()),k,0) for d in my_bus.ConnectedLoadsIDs], init=0)
        Q_ls = sum([get(get(grid.Loads[d].Q_shedding_tk,t,Dict()),k,0) for d in my_bus.ConnectedLoadsIDs], init=0)

        push!(grid.BusData_output,[bus_id, V, δ, Pg, Qg, Pd, Qd, P_ls, Q_ls, "AC"])
    end

    for bus_id in sort(collect(keys(grid.DCBuses)))
        my_bus = grid.DCBuses[bus_id]
        V = get(get(my_bus.V_magnitude_tk,t,Dict()),k,1)
        δ = 0
        Pg = sum([get(get(grid.DCGenerators[g].Pg_tk,t,Dict()),k,0) for g in my_bus.ConnectedGensIDs if grid.DCGenerators[g].GenType != :virtual], init=0)
        Qg = 0
        Pd = sum([grid.DCLoads[d].Pd_t[t] for d in my_bus.ConnectedLoadsIDs], init=0)
        Qd = 0
        P_ls = sum([get(get(grid.DCLoads[d].P_shedding_tk,t,Dict()),k,0) for d in my_bus.ConnectedLoadsIDs], init=0)
        Q_ls = 0
        push!(grid.BusData_output,[bus_id, V, δ, Pg, Qg, Pd, Qd, P_ls, Q_ls, "DC"])
    end

    # Populate "LineLoading" dataframe
    # BranchID = Int64[], FromBus = Int64[],ToBus = Int64[],PL_1 = Float64[],PL_2 = Float64[],PLoss = Float64[],QL_1 = Float64[],QL_2 = Float64[],QLoss = Float64[],Utilization = Float64[]
    for line_id in sort(collect(keys(grid.Branches)))
        my_branch = grid.Branches[line_id]
        FromBus = my_branch.Fr_bus_ID
        ToBus = my_branch.To_bus_ID
        
        PL_1 = get(get(my_branch.PowerFlow_ij_tk,t,Dict()),k,0)
        PL_2 = get(get(my_branch.PowerFlow_ji_tk,t,Dict()),k,0)
        PLoss = abs(PL_1 + PL_2)
        QL_1 = get(get(my_branch.ReactFlow_ij_tk,t,Dict()),k,0)
        QL_2 = get(get(my_branch.ReactFlow_ji_tk,t,Dict()),k,0)
        QLoss = abs(QL_1 + QL_2)
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating*grid.S_base), sqrt(PL_2^2+QL_2^2)/(my_branch.rating*grid.S_base)])
        if my_branch.BranchType == 0
            label = "AC-Branch"
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

        PL_1 = get(get(my_branch.PowerFlow_ij_tk,t,Dict()),k,0)
        PL_2 = get(get(my_branch.PowerFlow_ji_tk,t,Dict()),k,0)
        PLoss = abs(PL_1 + PL_2)
        QL_1 = 0
        QL_2 = 0
        QLoss = 0
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating*grid.S_base), sqrt(PL_2^2+QL_2^2)/(my_branch.rating*grid.S_base)])
        push!(grid.LineLoading,[line_id,FromBus,ToBus,PL_1,PL_2,PLoss,QL_1,QL_2,QLoss,Utilization*100, "DC-Branch"])
    end

    for line_id in sort(collect(keys(grid.DCLinks)))
        my_branch = grid.DCLinks[line_id]
        FromBus = my_branch.Fr_bus_ID
        ToBus = my_branch.To_bus_ID
        PL_1 = get(get(my_branch.PowerFlow_ij_tk,t,Dict()),k,0)
        PL_2 = get(get(my_branch.PowerFlow_ji_tk,t,Dict()),k,0)
        PLoss = abs(PL_1 + PL_2)
        QL_1 = 0
        QL_2 = 0
        QLoss = 0
        Utilization = maximum([sqrt(PL_1^2+QL_1^2)/(my_branch.rating), sqrt(PL_2^2+QL_2^2)/(my_branch.rating)])
        push!(grid.LineLoading,[line_id,FromBus,ToBus,PL_1,PL_2,PLoss,QL_1,QL_2,QLoss,Utilization*100, "DC-Link"])
    end

    # Flow through Converters
    for conv in sort(collect(keys(grid.Converters)))
        ConverterID = grid.Converters[conv].Conv_ID
        ac_gen_id = grid.Converters[conv].gen_ac_id
        dc_gen_id = grid.Converters[conv].gen_dc_id
        if grid.Converters[conv].type == :ACDC
            AC_Bus = grid.Generators[ac_gen_id].GenBus_ID
            DC_Bus = grid.DCGenerators[dc_gen_id].GenBus_ID
            P_ACDC = get(get(grid.Generators[ac_gen_id].Pg_tk,t,Dict()),k,0)
            P_DCAC = get(get(grid.DCGenerators[dc_gen_id].Pg_tk,t,Dict()),k,0)
        elseif grid.Converters[conv].type == :B2B
            AC_Bus = grid.Generators[ac_gen_id].GenBus_ID
            DC_Bus = grid.Generators[dc_gen_id].GenBus_ID
            P_ACDC = get(get(grid.Generators[ac_gen_id].Pg_tk,t,Dict()),k,0)
            P_DCAC = get(get(grid.Generators[dc_gen_id].Pg_tk,t,Dict()),k,0)
        end
        PLoss = abs(P_ACDC + P_DCAC)
        Utilization = maximum([abs(P_ACDC)/grid.Converters[conv].rate,abs(P_ACDC)/grid.Converters[conv].rate])
        type = string(grid.Converters[conv].type)
        push!(grid.Converter_flow, [ConverterID, AC_Bus, DC_Bus, P_ACDC,P_DCAC,PLoss,Utilization*100,type])
    end
end

function build_full_DOPF_model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    # Sanity check
    solver = SimulationSettings.MILP_solver
    NLP_flag = false
    if SimulationSettings.ac_grid_model == :AC || SimulationSettings.dc_grid_model == :FNL
        solver = SimulationSettings.NLP_solver
        NLP_flag = true
    end

    if NLP_flag
        if SimulationSettings.transmission_switching != [] || SimulationSettings.substation_switching != []
            error("Switching is not supported with Nonlinear problem settings.")
            return -1
        end
    end

    # Model initialization
    model = Model(solver)

    DOPF_variable_initialization!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_nodal_balance_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_ac_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_ac_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_ac_branch_dynamic!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_ac_branch_dynamic!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_substation_switching_ac_node!(model, grid, SimulationSettings, prerequisites_data)    
    DOPF_generator_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_load_shedding_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_voltage_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_angle_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data) 
    DOPF_converter_flow!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_converter_capacity!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_nodal_balance_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_dc_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_dc_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_dc_link!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_dc_link!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_voltage_limits_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_generator_limits_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_transition_constraints!(model, grid, SimulationSettings, prerequisites_data)
    if length(keys(prerequisites_data.fixed_schedules)) != 0
        DOPF_schedule_fixes!(model, grid, SimulationSettings, prerequisites_data)
    end
    if length(keys(prerequisites_data.fixed_topology)) != 0
        DOPF_substation_switching_ac_node_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_powerflow_ac_branch_dynamic_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_thermal_limits_ac_branch_dynamic_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_topology_fixes!(model, grid, SimulationSettings, prerequisites_data)
    end
    if length(keys(prerequisites_data.fixed_commitments)) != 0
        DOPF_commitment_fixes!(model, grid, SimulationSettings, prerequisites_data)
    end
    DOPF_objective_function!(model, grid, SimulationSettings, prerequisites_data)
 
    return model
end

function build_snapshot_DOPF_model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites; k=1)

    # Model initialization
    model = Model(solver)

    DOPF_variable_initialization!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_nodal_balance_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_ac_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_ac_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_ac_branch_dynamic!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_ac_branch_dynamic!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_substation_switching_ac_node!(model, grid, SimulationSettings, prerequisites_data)    
    DOPF_generator_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_load_shedding_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_voltage_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_angle_limits_ac_node!(model, grid, SimulationSettings, prerequisites_data) 
    DOPF_converter_flow!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_converter_capacity!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_nodal_balance_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_dc_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_dc_branch_static!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_powerflow_dc_link!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_thermal_limits_dc_link!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_voltage_limits_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    DOPF_generator_limits_dc_node!(model, grid, SimulationSettings, prerequisites_data)
    if length(keys(prerequisites_data.fixed_schedules)) != 0
        DOPF_schedule_fixes!(model, grid, SimulationSettings, prerequisites_data; k=k)
    end
    if length(keys(prerequisites_data.fixed_topology)) != 0
        DOPF_substation_switching_ac_node_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_powerflow_ac_branch_dynamic_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_thermal_limits_ac_branch_dynamic_fixed!(model, grid, SimulationSettings, prerequisites_data)
        DOPF_topology_fixes!(model, grid, SimulationSettings, prerequisites_data)
    end
    if length(keys(prerequisites_data.fixed_commitments)) != 0
        DOPF_commitment_fixes!(model, grid, SimulationSettings, prerequisites_data)
    end
    DOPF_objective_function!(model, grid, SimulationSettings, prerequisites_data)
 
    return model
end

function optimize_DOPF_model!(model::Model, grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook; update_grid=false, update_order_book=false)

    optimize!(model)

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

function calculate_opex_t(model::Model,grid::PowerGrid, prerequisites_data::DOPF_Prerequisites, t; include_load_shedding=true, include_commitment_cost=true)
    if t ∉ prerequisites_data.time_horizon
        error("Selected time instance ($t) is out of bounds")
    else
        GenBids = prerequisites_data.Order_Book.Gen_bids
        LoadBids = prerequisites_data.Order_Book.Load_bids
        
        val_p_gen_ac = sum([JuMP.value.(model[:p_gen_ac])[g,1,t]*GenBids[g]["price"][t][1] for g in prerequisites_data.ac_gen_ids], init=0)
        val_commitment = sum([grid.Generators[g].C0*JuMP.value.(model[:u_gt])[g,t]+JuMP.value.(model[:α_gt])[g,t]*grid.Generators[g].start_up_cost+JuMP.value.(model[:β_gt])[g,t]*grid.Generators[g].shut_down_cost for g in prerequisites_data.commitable_gen_ids], init=0)
        
        val_p_gen_dc = sum([JuMP.value.(model[:p_gen_dc])[g,1,t]*grid.DCGenerators[g].C1 for g in prerequisites_data.dc_gen_ids], init=0)
        val_p_load_shedding = sum([JuMP.value.(model[:p_ls_ac])[d,k,t]*LoadBids[d]["price"][t] for d  in prerequisites_data.ac_load_shedding_ids, k in prerequisites_data.k_t[t]], init=0)
        
        total = val_p_gen_ac + val_p_gen_dc

        if include_commitment_cost
            total += val_commitment
        end

        if include_load_shedding
            total += val_p_load_shedding
        end

        return total
    end
end

function clear_market!(market::PowerMarket, grid::PowerGrid, order_book::OrderBook, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites; update_grid=false)

    if order_book.Clearing == :UC # unit commitment -> grid-agnostic
        
        for t in keys(prerequisites_data.k_t)
            prerequisites_data.k_t[t] = [1]
        end

        model, solved_flag = UC_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        
        if solved_flag

            dual_model, solved_flag_relaxed = Fixed_UC_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
            
            if ! solved_flag_relaxed
                error("Primal relaxed problem is infeasibile!!")
            end

            if JuMP.has_duals(dual_model)
                μ_t = Dict()
                opex = Dict()
                for t in prerequisites_data.time_horizon
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => calculate_opex_t(dual_model,grid,prerequisites_data,t))
                end

                prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
                prerequisites_data.Order_Book.Opex = deepcopy(opex)
                prerequisites_data.time_horizon = simulation_settings.time_horizon
                settle_order_book!(grid, prerequisites_data)
                update_original_order_book!(order_book, prerequisites_data)
                compensate!(market, order_book)
                
            else
                error("Relaxed problem has no dual solution, which is strange if you ask me!")
            end

        end

        return solved_flag

    elseif order_book.Clearing == :NCUC # network constrained unit commitment -> grid-aware
        for t in keys(prerequisites_data.k_t)
            prerequisites_data.k_t[t] = [1]
        end
        model, solved_flag = NCUC_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        
        if solved_flag

            dual_model, solved_flag_relaxed = Fixed_SCUC_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
            
            if ! solved_flag_relaxed
                error("Primal relaxed problem is infeasibile!!")
            end

            if JuMP.has_duals(dual_model)
                μ_t = Dict()
                opex = Dict()
                lmp_t = Dict()
                for t in prerequisites_data.time_horizon
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => calculate_opex_t(dual_model,grid,prerequisites_data,t))
                    lmp_i = Dict()
                    for i in prerequisites_data.ac_node_ids
                        push!(lmp_i, i => sum(dual(dual_model[:ac_active_nodal_balance][i, k, t]) for k in prerequisites_data.k_t[t]))
                    end
                    push!(lmp_t, t => lmp_i)
                end

                prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
                prerequisites_data.Order_Book.Opex = deepcopy(opex)
                prerequisites_data.Order_Book.λ_i_t = deepcopy(lmp_t)
                prerequisites_data.time_horizon = simulation_settings.time_horizon
                settle_order_book!(grid, prerequisites_data)
                update_original_order_book!(order_book, prerequisites_data)
                compensate!(market, order_book)
                for gen_duplet_id in keys(grid.Generator_Duplets)
                    merge_element!(grid, :Generator, gen_duplet_id)
                end
            else
                error("Relaxed problem has no dual solution, which is strange if you ask me!")
            end

        end

        return solved_flag

    elseif order_book.Clearing == :SCUC # security constrained unit commitment -> grid-aware
        
        model, solved_flag = SCUC_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        
        if solved_flag

            dual_model, solved_flag_relaxed = Fixed_SCUC_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
            
            if ! solved_flag_relaxed
                error("Primal relaxed problem is infeasibile!!")
            end

            if JuMP.has_duals(dual_model)
                μ_t = Dict()
                opex = Dict()
                lmp_t = Dict()
                for t in prerequisites_data.time_horizon
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => calculate_opex_t(dual_model,grid,prerequisites_data,t))
                    lmp_i = Dict()
                    for i in prerequisites_data.ac_node_ids
                        push!(lmp_i, i => sum(dual(dual_model[:ac_active_nodal_balance][i, k, t]) for k in prerequisites_data.k_t[t]))
                    end
                    push!(lmp_t, t => lmp_i)
                end

                prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
                prerequisites_data.Order_Book.Opex = deepcopy(opex)
                prerequisites_data.Order_Book.λ_i_t = deepcopy(lmp_t)
                prerequisites_data.time_horizon = simulation_settings.time_horizon
                settle_order_book!(grid, prerequisites_data)
                update_original_order_book!(order_book, prerequisites_data)
                compensate!(market, order_book)
                for gen_duplet_id in keys(grid.Generator_Duplets)
                    merge_element!(grid, :Generator, gen_duplet_id)
                end
            else
                error("Relaxed problem has no dual solution, which is strange if you ask me!")
            end

        end

        return solved_flag
    elseif order_book.Clearing == :ED # economic dispatch -> grid-agnostic
        for t in keys(prerequisites_data.k_t)
            prerequisites_data.k_t[t] = [1]
        end

        model, solved_flag = ED_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        
        if solved_flag

            dual_model, solved_flag_relaxed = Fixed_UC_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
            
            if ! solved_flag_relaxed
                error("Primal relaxed problem is infeasibile!!")
            end

            if JuMP.has_duals(dual_model)
                μ_t = Dict()
                opex = Dict()
                for t in prerequisites_data.time_horizon
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => calculate_opex_t(dual_model,grid,prerequisites_data,t))
                end

                prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
                prerequisites_data.Order_Book.Opex = deepcopy(opex)
                prerequisites_data.time_horizon = simulation_settings.time_horizon
                settle_order_book!(grid, prerequisites_data)
                update_original_order_book!(order_book, prerequisites_data)
                compensate!(market, order_book)
                
            else
                error("Relaxed problem has no dual solution, which is strange if you ask me!")
            end

        end

        return solved_flag
    elseif order_book.Clearing == :OPF # optimal power flow -> grid-aware
        
        for t in keys(prerequisites_data.k_t)
            prerequisites_data.k_t[t] = [1]
        end
        solved_flag = true
        μ_t = Dict()
        opex = Dict()
        lmp_t = Dict()

        for t in simulation_settings.time_horizon
            prerequisites_data.time_horizon = [t]

            model, solved_flag = OPF_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        

            if solved_flag

                if JuMP.has_duals(model)
                    lmp_line = Dict()
                    for line_id in keys(grid.Branches)
                        if line_id in prerequisites_data.relaxed_capacity_lines
                            push!(lmp_line, line_id => 0)
                        else
                            a = sum(dual(model[:ac_active_flow_limits][line_id, 1, k, t]) for k in prerequisites_data.k_t[t])
                            b = sum(dual(model[:ac_active_flow_limits][line_id, 2, k, t]) for k in prerequisites_data.k_t[t])
                            push!(lmp_line, line_id => maximum([abs(a),abs(b)])) 
                        end
                    end

                    lmp_node = Dict()
                    for node_id in keys(grid.Buses)
                        push!(lmp_node, node_id => sum(dual(model[:ac_active_nodal_balance][node_id, k, t]) for k in prerequisites_data.k_t[t]))
                    end
                    push!(grid.Line_Duals, t => lmp_line)
                    push!(grid.Bus_Duals, t => lmp_node)
                end

                objective_value = JuMP.objective_value(model)
                dual_model, solved_flag_relaxed = Fixed_SCOPF_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
                
                if ! solved_flag_relaxed
                    error("Primal relaxed problem is infeasibile!!")
                end
    
                if JuMP.has_duals(dual_model)
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => objective_value)
                    lmp_i = Dict()
                    for i in prerequisites_data.ac_node_ids
                        push!(lmp_i, i => sum(dual(dual_model[:ac_active_nodal_balance][i, k, t]) for k in prerequisites_data.k_t[t]))
                    end
                    push!(lmp_t, t => lmp_i)
                else
                    error("Relaxed problem has no dual solution, which is strange if you ask me!")
                end
            else
                push!(μ_t, t => 0)
                push!(opex, t => -1)
                lmp_i = Dict()
                for i in prerequisites_data.ac_node_ids
                    push!(lmp_i, i => 0)
                end
                push!(lmp_t, t => lmp_i)
                println("================================")
                println("Couldn't solve time instance: $t")
                println("================================")
            end
        end
        prerequisites_data.time_horizon = simulation_settings.time_horizon

        prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
        prerequisites_data.Order_Book.Opex = deepcopy(opex)
        prerequisites_data.Order_Book.λ_i_t = deepcopy(lmp_t)

        settle_order_book!(grid, prerequisites_data)
        update_original_order_book!(order_book, prerequisites_data)
        compensate!(market, order_book)
        for gen_duplet_id in keys(grid.Generator_Duplets)
            merge_element!(grid, :Generator, gen_duplet_id)
        end
        return solved_flag
    elseif order_book.Clearing == :SCOPF # security constrained optimal power flow -> grid-aware

        solved_flag = true
        μ_t = Dict()
        opex = Dict()
        lmp_t = Dict()

        for t in simulation_settings.time_horizon
            prerequisites_data.time_horizon = [t]

            model, solved_flag = SCOPF_Model!(grid, simulation_settings, prerequisites_data, order_book ; update_grid=update_grid)
        
            if solved_flag

                if JuMP.has_duals(model)
                    lmp_line = Dict()
                    for line_id in keys(grid.Branches)
                        if line_id in prerequisites_data.relaxed_capacity_lines
                            push!(lmp_line, line_id => 0)
                        else
                            push!(lmp_line, line_id => sum(dual(model[:ac_active_flow_limits][line_id, 1, k, t]) for k in prerequisites_data.k_t[t]))
                        end
                    end

                    lmp_node = Dict()
                    for node_id in keys(grid.Buses)
                        a = sum(dual(model[:ac_active_flow_limits][line_id, 1, k, t]) for k in prerequisites_data.k_t[t])
                        b = sum(dual(model[:ac_active_flow_limits][line_id, 2, k, t]) for k in prerequisites_data.k_t[t])
                        push!(lmp_line, line_id => maximum([abs(a),abs(b)])) 
                    end
                    push!(grid.Line_Duals, t => lmp_line)
                    push!(grid.Bus_Duals, t => lmp_node)
                end

                objective_value = JuMP.objective_value(model)
                dual_model, solved_flag_relaxed = Fixed_SCOPF_Model!(model, grid, simulation_settings, prerequisites_data, order_book)
                
                if ! solved_flag_relaxed
                    error("Primal relaxed problem is infeasibile!!")
                end
    
                if JuMP.has_duals(dual_model)
                    push!(μ_t, t => sum([JuMP.dual(dual_model[:single_node_balance][k,t]) for k in prerequisites_data.k_t[t]], init=0))
                    push!(opex, t => objective_value)
                    lmp_i = Dict()
                    for i in prerequisites_data.ac_node_ids
                        push!(lmp_i, i => sum(dual(dual_model[:ac_active_nodal_balance][i, k, t]) for k in prerequisites_data.k_t[t]))
                    end
                    push!(lmp_t, t => lmp_i)
                else
                    error("Relaxed problem has no dual solution, which is strange if you ask me!")
                end
            else
                push!(μ_t, t => 0)
                push!(opex, t => -1)
                lmp_i = Dict()
                for i in prerequisites_data.ac_node_ids
                    push!(lmp_i, i => 0)
                end
                push!(lmp_t, t => lmp_i)
                println("================================")
                println("Couldn't solve time instance: $t")
                println("================================")
            end
        end
        prerequisites_data.time_horizon = simulation_settings.time_horizon

        prerequisites_data.Order_Book.μ_t = deepcopy(μ_t)
        prerequisites_data.Order_Book.Opex = deepcopy(opex)
        prerequisites_data.Order_Book.λ_i_t = deepcopy(lmp_t)

        settle_order_book!(grid, prerequisites_data)
        update_original_order_book!(order_book, prerequisites_data)
        compensate!(market, order_book)
        for gen_duplet_id in keys(grid.Generator_Duplets)
            merge_element!(grid, :Generator, gen_duplet_id)
        end
        return solved_flag
    end
end

function settle_order_book!(grid::PowerGrid, prerequisites_data::DOPF_Prerequisites)
    order_book = prerequisites_data.Order_Book
    if order_book.Pricing == :Uniform
        # generators ramping up (+ve output) get payed μ_t while generators ramping down (-ve output) get payed -C1 (pay their marginal costs)
        for gen_id in keys(order_book.Schedule["gen"])
            for t in prerequisites_data.time_horizon
                if order_book.Schedule["gen"][gen_id][t] ≥ 0
                    push!(order_book.Revenue[gen_id], t => order_book.Schedule["gen"][gen_id][t]*order_book.μ_t[t])
                else
                    push!(order_book.Revenue[gen_id], t => order_book.Schedule["gen"][gen_id][t]*grid.Generators[gen_id].C1)
                end
            end
        end

        for load_id in keys(order_book.Schedule["load"])
            for t in prerequisites_data.time_horizon
                load_payment = (order_book.Load_bids[load_id]["qty"][t]+order_book.Schedule["load"][load_id][t])*order_book.μ_t[t]
                load_revenue = order_book.Schedule["load"][load_id][t]*order_book.Load_bids[load_id]["price"][t]
                push!(order_book.Payment[load_id], t => load_payment-load_revenue)
            end
        end
    elseif order_book.Pricing == :PAB
        # All generators get payed what they ask for
        for gen_id in keys(order_book.Schedule["gen"])
            for t in prerequisites_data.time_horizon
                push!(order_book.Revenue[gen_id], t => order_book.Schedule["gen"][gen_id][t]*order_book.Gen_bids[gen_id]["price"][t])
            end
        end
        for t in prerequisites_data.time_horizon
            ΣLS = sum([order_book.Schedule["load"][d][t] for d in keys(order_book.Schedule["load"])])
            ΣPₗ = sum([order_book.Load_bids[d]["qty"][t] for d in keys(order_book.Schedule["load"])]) + ΣLS
            Σπ = sum([order_book.Revenue[g][t] for g in keys(order_book.Schedule["gen"])])
            AVG_cost = Σπ/ΣPₗ
            for load_id in keys(order_book.Schedule["load"])
                load_payment = (order_book.Load_bids[load_id]["qty"][t]+order_book.Schedule["load"][load_id][t])*AVG_cost
                load_revenue = order_book.Schedule["load"][load_id][t]*order_book.Load_bids[load_id]["price"][t]
                push!(order_book.Payment[load_id], t => load_payment-load_revenue)
            end
        end
    elseif order_book.Pricing == :LMP
        # generators ramping up (+ve output) get payed λ_i_t while generators ramping down (-ve output) get payed -C1 (pay their marginal costs)
    end
    prerequisites_data.Order_Book = order_book
end

function update_original_order_book!(original_order_book::OrderBook, prerequisites_data::DOPF_Prerequisites)
    for t in prerequisites_data.time_horizon
        for g in prerequisites_data.Order_Book.Gen_ids
            g_root = prerequisites_data.ac_gen_id_to_gen_root[g]
            original_order_book.Schedule["gen"][g_root][t] += prerequisites_data.Order_Book.Schedule["gen"][g][t]
            original_order_book.Revenue[g_root][t] += prerequisites_data.Order_Book.Revenue[g][t]
        end

        for d in prerequisites_data.ac_load_shedding_ids
            original_order_book.Schedule["load"][d][t] += prerequisites_data.Order_Book.Schedule["load"][d][t]
            original_order_book.Payment[d][t] += prerequisites_data.Order_Book.Payment[d][t]
        end
    end
    original_order_book.Opex  = deepcopy(prerequisites_data.Order_Book.Opex)
    original_order_book.μ_t   = deepcopy(prerequisites_data.Order_Book.μ_t)
    original_order_book.λ_i_t = deepcopy(prerequisites_data.Order_Book.λ_i_t)
end

function DOPF_single_node_model!(grid::PowerGrid, simulation_settings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    # Model initialization
    if simulation_settings.ac_grid_model == :Bθ
        solver = simulation_settings.MILP_solver
    else
        solver = simulation_settings.NLP_solver
    end
    model = Model(solver)
    DOPF_single_node_variable_initialization!(model, grid, simulation_settings, prerequisites_data)
    DOPF_single_node_balance!(model, grid, simulation_settings, prerequisites_data)
    DOPF_generator_limits_ac_node!(model, grid, simulation_settings, prerequisites_data)
    DOPF_load_shedding_limits_ac_node!(model, grid, simulation_settings, prerequisites_data)
    DOPF_generator_limits_dc_node!(model, grid, simulation_settings, prerequisites_data)
    DOPF_transition_constraints!(model, grid, simulation_settings, prerequisites_data)
    DOPF_schedule_fixes!(model, grid, simulation_settings, prerequisites_data)
    DOPF_topology_fixes!(model, grid, simulation_settings, prerequisites_data)
    DOPF_commitment_fixes!(model, grid, simulation_settings, prerequisites_data)
    DOPF_objective_function!(model, grid, simulation_settings, prerequisites_data)
    return model
end

function UC_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    model = DOPF_single_node_model!(grid, SimulationSettings, prerequisites_data)
    solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)

    return model, solution_status
end

function NCUC_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    model = build_full_DOPF_model!(grid, SimulationSettings, prerequisites_data)
    solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)
    return model, solution_status
end

function SCUC_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    if SimulationSettings.Meta_solver == :none
        model = build_full_DOPF_model!(grid, SimulationSettings, prerequisites_data)
        solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)
        return model, solution_status
    elseif SimulationSettings.Meta_solver == :CCG
        # Do something else
    end
end

function ED_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    model = DOPF_single_node_model!(grid, SimulationSettings, prerequisites_data)
    solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)

    return model, solution_status
end

function OPF_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    if SimulationSettings.Meta_solver == :none
        model = build_full_DOPF_model!(grid, SimulationSettings, prerequisites_data)
        solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)
        return model, solution_status
    elseif SimulationSettings.Meta_solver == :CCG

    end
end

function SCOPF_Model!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false)
    if SimulationSettings.Meta_solver == :none
        model = build_full_DOPF_model!(grid, SimulationSettings, prerequisites_data)
        solution_status = optimize_DOPF_model!(model, grid, SimulationSettings, prerequisites_data, order_book; update_grid=update_grid, update_order_book=true)
        return model, solution_status
    elseif SimulationSettings.Meta_solver == :CCG

    end
end

