using Parameters

@with_kw mutable struct OrderBook
    Market = "day-ahead" # could be any market
    Clearing = [] # supports :UC and :ED for zonal schemes and :NCUC - :SCUC - :OPF - :SCOPF for nodal schemes
    Pricing = [] # supports :PAB (pay as bid) - :LMP (Locational marginal pricing) - :Uniform 
    horizon = 1:24
    Gen_ids = []
    Load_ids = []
    Gen_bids = Dict() # dictionary of generator bids --- [gen_id -> qty/price -> t -> (val_up,val_down)]
    Load_bids = Dict() # dictionary of load bids --- [load_id -> qty/price -> t -> val] --- price and qty here mean amounts for shedding
    Schedule = Dict("gen"=>Dict(), "load"=>Dict()) # dictionary of schedules after market clearing ["gen"/"load" -> gen_id/load_id -> t -> P_gt/P_dt]
    Revenue = Dict() # dictionary with each generator's revenue
    Payment = Dict() # dictionary of payments by each load
    Opex = Dict() # hourly operating cost -> after day-ahead market clearing
    μ_t = Dict() # market clearing price for each time step
    λ_i_t = Dict() # LMPs for each node and each time step
    Grid_Options = Dict() # contains information on substations, relaxed physics, or relaxed technical constraints, if applicable
    Grid_Results = Dict() # contains information about grid results after simulation
end

@with_kw mutable struct PowerMarket

    Name = "MyMarket" # can be anything, even the date of clearing
    horizon = 1:24 # time steps the market is solved for at each day
    Gen_ids = [] # real generators
    Load_ids = [] # real loads
    Grid = [] # grid is required in case of nodal markets where power flows need to be calculated
    Order_Books = [] # List of order books for different markets -> 2*D
    Schedule_list = [] # list of schedules for different days -> D
    Revenue_list = []
    Payment_list = []
    Opex_list = []
    μ_t_list = []
    λ_i_t_list = []

    Schedule = Dict("gen"=>Dict(), "load"=>Dict()) # dictionary of schedules after market clearing ["gen"/"load" -> gen_id/load_id -> t -> P_gt/P_dt]
    Revenue = Dict() # dictionary with each generator's revenue
    Payment = Dict() # dictionary of payments by each load
    Opex = Dict() # hourly operating cost -> after day-ahead market clearing
    μ_t = Dict() # market clearing price for each time step
    λ_i_t = Dict() # LMPs for each node and each time step
end

function create_market(Gen_ids, Load_ids; name="MyMarket", horizon=1:24)
    my_market = PowerMarket(Gen_ids=Gen_ids, Load_ids=Load_ids, Name=name, horizon=horizon) 
    setup_market!(my_market)
    return my_market  
end

function attach_grid2market!(market ::PowerMarket, grid ::PowerGrid)
    market.Grid = deepcopy(grid)
    return nothing
end

function setup_market!(market ::PowerMarket)
    schedule_t = Dict()
    revenue_t = Dict()
    for t in market.horizon
        push!(schedule_t, t => 0)
        push!(revenue_t, t => 0)
        push!(market.Opex, t => 0)
    end
    
    for gen_id in keys(market.Gen_ids)
        push!(market.Schedule["gen"], gen_id => deepcopy(schedule_t))
        push!(market.Revenue, gen_id => deepcopy(revenue_t))
    end
    for load_id in keys(market.Load_ids)
        push!(market.Schedule["load"], load_id => deepcopy(schedule_t))
        push!(market.Payment, load_id => deepcopy(revenue_t))
    end

    return nothing
end

function initialize_session!(order_book ::OrderBook)
    """
    This function is called just after the order book is created
    """
    schedule_t = Dict()
    revenue_t = Dict()
    for t in order_book.horizon
        push!(schedule_t, t => 0)
        push!(revenue_t, t => 0)
    end

    for gen_id in keys(order_book.Gen_ids)
        push!(order_book.Schedule["gen"], gen_id => deepcopy(schedule_t))
        push!(order_book.Revenue, gen_id => deepcopy(revenue_t))
    end
    for load_id in keys(order_book.Load_ids)
        push!(order_book.Schedule["load"], load_id => deepcopy(schedule_t))
        push!(order_book.Payment, load_id => deepcopy(revenue_t))
    end
end

function compensate!(market ::PowerMarket, order_book ::OrderBook)

    for t in collect(market.horizon)
        market.Opex[t] += order_book.Opex[t]
    end

    for gen_id in keys(order_book.Schedule["gen"])
        for t in collect(market.horizon)
            market.Revenue[gen_id][t] += order_book.Revenue[gen_id][t]
            market.Schedule["gen"][gen_id][t] += order_book.Schedule["gen"][gen_id][t]
        end
    end

    for load_id in keys(order_book.Schedule["load"])
        for t in collect(market.horizon)
            market.Payment[load_id][t] += order_book.Schedule["load"][load_id][t]*order_book.μ_t[t]
            market.Schedule["load"][load_id][t] += order_book.Schedule["load"][load_id][t]
        end
    end

    push!(market.Order_Books, order_book)
end

function submit_gen_bid!(order_book ::OrderBook, gen_id, price, qty)
    # price is a dictionary for offering price at each time step
    # qty is a dictionary for offering quantity at each time step
    bid = Dict("price"=>price, "qty"=>qty)
    push!(order_book.Gen_bids, gen_id => bid)    
end

function submit_load_bid!(order_book ::OrderBook, load_id, price, qty)
    # qty is a dictionary of the required load demand at each time step
    bid = Dict("price"=>price, "qty"=>qty)
    push!(order_book.Load_bids, load_id => bid)
end

function create_bids_from_grid!(market ::PowerMarket, order_book ::OrderBook)
    
    if market.Grid == []
        error("No grid attached! Use `attach_grid2market!(market ::PowerMarket, grid ::PowerGrid)`")
    end

    for gen_id in market.Gen_ids
        qty_gen = Dict()
        price_gen = Dict()
        for t in market.horizon
            push!(qty_gen, t => [market.Grid.Generators[gen_id].Pg_max, market.Grid.Generators[gen_id].Pg_min])
            push!(price_gen, t => [market.Grid.Generators[gen_id].C1])
        end
        submit_gen_bid!(order_book, gen_id, price_gen, qty_gen)
    end

    for load_id in market.Load_ids
        qty_load = Dict()
        price_load = Dict()
        
        t_total = length(market.Grid.Loads[load_id].Pd_t)
        for t in market.horizon
            if t > t_total
                push!(qty_load, t => market.Grid.Loads[load_id].Pd_t[end])
            else
                push!(qty_load, t => market.Grid.Loads[load_id].Pd_t[t])
            end
            push!(price_load, t => market.Grid.Loads[load_id].Shedding_Cost)
        end
        submit_load_bid!(order_book, load_id, price_load, qty_load)
    end
end

function create_auto_base_bids!(market::PowerMarket, order_book::OrderBook; new_day=true)
    if new_day
        push!(market.Schedule_list, deepcopy(market.Schedule))
        push!(market.Revenue_list, deepcopy(market.Revenue))
        push!(market.Payment_list, deepcopy(market.Payment))
        push!(market.Opex_list, deepcopy(market.Opex))
        push!(market.μ_t_list, deepcopy(market.μ_t))
        push!(market.λ_i_t_list, deepcopy(market.λ_i_t))

        setup_market!(market)
        create_bids_from_grid!(market, order_book)
    else
        create_auto_RD_bids!(market, order_book)
    end
end

function create_auto_RD_bids!(market ::PowerMarket, order_book ::OrderBook)
    for gen_id in market.Gen_ids
        qty_gen = Dict() # (qty_up, qty_down)
        price_gen = Dict()
        for t in market.horizon
            push!(qty_gen, t => [market.Grid.Generators[gen_id].Pg_max - market.Schedule["gen"][gen_id][t], -market.Schedule["gen"][gen_id][t] + market.Grid.Generators[gen_id].Pg_min])
            push!(price_gen, t => [market.Grid.Generators[gen_id].C1*2, market.Grid.Generators[gen_id].C1])
        end
        submit_gen_bid!(order_book, gen_id, price_gen, qty_gen)
    end

    for load_id in market.Load_ids
        qty_load = Dict()
        price_load = Dict()
        
        t_total = length(market.Grid.Loads[load_id].Pd_t)
        for t in market.horizon
            if t > t_total
                push!(qty_load, t => market.Grid.Loads[load_id].Pd)
            else
                push!(qty_load, t => market.Grid.Loads[load_id].Pd_t[t])
            end
            push!(price_load, t => market.Grid.Loads[load_id].Shedding_Cost*10)
        end
        submit_load_bid!(order_book, load_id, price_load, qty_load)
    end
end

function add_sided_generators!(market ::PowerMarket)
    for gen_id in deepcopy(keys(market.Grid.Generators))
        split_element!(market.Grid, :Generator, gen_id)
    end
end

# function ED!(schedule::Dict, order_book::OrderBook, t)
#     Gen_set = order_book.Gen_ids
#     Load_set = order_book.Load_ids

#     model = JuMP.Model(Gurobi.Optimizer)
    
#     # 1. Variables
#     JuMP.@variable(model, 0 ≤ p[g in Gen_set] ≤ get(order_book.Gen_bids[g]["qty"], t, 0))
#     # 2. Constraints
#     JuMP.@constraint(model,Nodal_balance,
#         sum(p[g] for g in Gen_set) + sum(schedule["gen"][g][t] for g in Gen_set) == sum(get(order_book.Load_bids[d]["qty"], t, 0) for d in Load_set) + sum(schedule["load"][d][t] for d in Load_set))

#     # 3. Objective
#     JuMP.@objective(model, Min, sum(get(order_book.Gen_bids[g]["price"], t, 0)*p[g] for g in Gen_set))
#     JuMP.optimize!(model)

#     MCP = JuMP.dual(model[:Nodal_balance])

#     Pg = JuMP.value.(model[:p])

#     for gen_id in Gen_set
#         push!(order_book.Schedule["gen"][gen_id], t => Pg[gen_id])
#     end

#     for load_id in Load_set
#         push!(order_book.Schedule["load"][load_id], t => get(order_book.Load_bids[load_id]["qty"], t, 0))
#     end
#     opex = JuMP.objective_value(model)
#     return MCP, opex
# end

# function SCOPF!(grid::PowerGrid, schedule ::Dict, order_book ::OrderBook, t)
#     # TODO: Continue work here 
#     Substations = get(order_book.Grid_Options, "substations", [])
#     Switchable_Branches = get(order_book.Grid_Options, "switchable_branches", [])
#     Transmission_Switching = get(order_book.Grid_Options, "transmission_switching", [])
#     Substation_Switching = get(order_book.Grid_Options, "substation_switching", Dict("reconf" => [], "splitting" => []))
#     no_load_shedding_pre_contingency = get(order_book.Grid_Options, "no_load_shedding_pre", true)
#     max_ts = Inf
#     max_reconf = Inf
#     max_splitting = Inf
    
#     # transmission_switching = 1
#     # if Transmission_Switching == []
#     #     transmission_switching = 0
#     # end

#     # substation_switching = 1
#     # if Substation_Switching == []
#     #     substation_switching = 0
#     # end

#     prerequisites_data = compile_simulation_prerequisites_SC!(grid, :Transmission_Contingencies, Switchable_Branches, Substations, [], Dict())

#     # activating active control of loads for load shedding
#     for load_id in keys(grid.Loads)
#         grid.Loads[load_id].ActiveControl = true
#     end

#     simulation_settings = SCOPF_SimulationSettings(ac_grid_model=:DCOPF,
#         transmission_switching=Transmission_Switching,
#         substation_switching=Substation_Switching,
#         max_transmission_switching=Inf, max_substation_reconf=Inf, max_busbar_splitting=Inf,
#         redispatch=false,
#         contingency_type=:Transmission_Contingencies,
#         NLP_solver=Ipopt.Optimizer, MILP_solver=Gurobi.Optimizer, Meta_solver=nothing)

#     # Build the model here
#     model = Model(Gurobi.Optimizer)

#     JuMP.@variable(model, δ[i in prerequisites_data.nodes_set, k in prerequisites_data.k])
#     JuMP.@variable(model, pij[i in prerequisites_data.nodes_set, j in prerequisites_data.nodes_set, k in prerequisites_data.k; Set([i, j]) in prerequisites_data.branch_nodes])
#     JuMP.@variable(model, p[g in prerequisites_data.gen_set, k in prerequisites_data.k])
#     JuMP.@variable(model, p_ls[d in prerequisites_data.load_set, k in prerequisites_data.k; grid.Loads[d].ActiveControl == true])

#     # Switching variables
#     # if prerequisites_data.switched_transmission_set != []
#     JuMP.@variable(model, z[l in prerequisites_data.switched_transmission_set, k in prerequisites_data.k], Bin)
#     # end

#     # Sanity check
#     if simulation_settings.substation_switching["reconf"] == [:post] && simulation_settings.substation_switching["splitting"] == [:pre]
#         error("Wrong problem setting!")
#         return -1
#     elseif simulation_settings.substation_switching["reconf"] == [:pre,:post] && length(simulation_settings.substation_switching["splitting"]) == 1
#         error("Wrong problem setting!")
#         return -1
#     elseif simulation_settings.substation_switching["splitting"] == [:pre,:post] && simulation_settings.substation_switching["reconf"] == [:post]
#         error("Wrong problem setting!")
#         return -1
#     end

#     # if prerequisites_data.reconf_set != []
#     JuMP.@variable(model, z_l[l in prerequisites_data.reconf_set, k in prerequisites_data.k], Bin)
#     # end

#     # if prerequisites_data.coupler_set != []
#     JuMP.@variable(model, z_c[c in prerequisites_data.coupler_set, k in prerequisites_data.k], Bin)
#     # end

#     #======================= Constriants =======================#

#     up_gen_set = []
#     down_gen_set = []
#     gen_id_dict = Dict()

#     for duplet_id in keys(grid.Generator_Duplets)
#         push!(up_gen_set, grid.Generator_Duplets[duplet_id][1])
#         push!(down_gen_set, grid.Generator_Duplets[duplet_id][2])
#         push!(gen_id_dict, grid.Generator_Duplets[duplet_id][1] => grid.Generator_Duplets[duplet_id][1])
#         push!(gen_id_dict, grid.Generator_Duplets[duplet_id][2] => grid.Generator_Duplets[duplet_id][1])
#     end

#     Nodes_set = prerequisites_data.nodes_set
#     Branch_nodes = prerequisites_data.branch_nodes
#     load_set = prerequisites_data.load_set
#     Gen_set = prerequisites_data.gen_set

#     JuMP.@constraint(model, Pnodal[i in Nodes_set,k in prerequisites_data.k],
#             sum(model[:pij][i,j,k] for j = Nodes_set if Set([i,j]) in Branch_nodes) == sum(model[:p][g,k] for g in Gen_set if grid.Generators[g].GenBus_ID == i) 
#                 + sum(schedule["gen"][g][t] for g in order_book.Gen_ids if grid.Generators[g].GenBus_ID == i) 
#                 - sum(schedule["load"][d][t] for d in prerequisites_data.load_set if grid.Loads[d].LoadBus_ID == i)
#                 + sum(model[:p_ls][l,k] for l in prerequisites_data.load_set if grid.Loads[l].LoadBus_ID == i))
    
#     Nodes_set = prerequisites_data.nodes_set
#     Branch_nodes = prerequisites_data.branch_nodes
#     unswitched_Transmission_nodes = prerequisites_data.unswitched_Transmission_nodes
#     branch_dictionary = prerequisites_data.branch_dictionary
#     Sbase = prerequisites_data.base_MVA
#     B = prerequisites_data.B_matrix

#     # ACTIVE POWER THROUGH LINE i-j
#     JuMP.@constraint(model,pl[i in Nodes_set,j in Nodes_set, k in prerequisites_data.k; Set([i,j]) in unswitched_Transmission_nodes],
#         model[:pij][i,j,k] == Sbase*(B[i,j])*(model[:δ][i,k]-model[:δ][j,k])*prerequisites_data.a_l[i,j,k])

#     # OPPOSITE FLOW CONSISTENCY
#     JuMP.@constraint(model, pl_consistency[i in Nodes_set,j in Nodes_set, k in prerequisites_data.k; Set([i,j]) in Branch_nodes], 
#         model[:pij][i,j,k] == -model[:pij][j,i,k])

#     Nodes_set = prerequisites_data.nodes_set
#     Sbase = prerequisites_data.base_MVA
#     branch_dictionary = prerequisites_data.branch_dictionary
#     unswitched_Transmission_nodes = prerequisites_data.unswitched_Transmission_nodes

#     # LINE CAPACITY
#     JuMP.@constraint(model,pl_rate[i in Nodes_set,j in Nodes_set, k in prerequisites_data.k; Set([i,j]) in unswitched_Transmission_nodes],
#         -Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*prerequisites_data.a_l[i,j,k] ≤ model[:pij][i,j,k] ≤ Sbase*grid.Branches[branch_dictionary[Set([i,j])]].rating*prerequisites_data.a_l[i,j,k])


#     Sbase = prerequisites_data.base_MVA
#     Nodes_set = prerequisites_data.nodes_set
#     switched_Transmission_nodes = prerequisites_data.switched_Transmission_nodes
#     branch_dictionary = prerequisites_data.branch_dictionary
#     # LINE CAPACITY
#     JuMP.@constraint(model, pl_rate_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in switched_Transmission_nodes],
#         -Sbase * grid.Branches[branch_dictionary[Set([i, j])]].rating * model[:z][branch_dictionary[Set([i, j])],k] * prerequisites_data.a_l[i,j,k] ≤ model[:pij][i, j, k])

#     JuMP.@constraint(model, pl_rate_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in switched_Transmission_nodes],
#         model[:pij][i, j, k] ≤ Sbase * grid.Branches[branch_dictionary[Set([i, j])]].rating * model[:z][branch_dictionary[Set([i, j])],k] * prerequisites_data.a_l[i,j,k])

#     switched_Transmission_nodes = prerequisites_data.switched_Transmission_nodes
#     Nodes_set = prerequisites_data.nodes_set
#     branch_dictionary = prerequisites_data.branch_dictionary
#     B = prerequisites_data.B_matrix
#     Sbase = prerequisites_data.base_MVA
#     M = prerequisites_data.base_MVA * 100
#     JuMP.@constraint(model, pl_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in switched_Transmission_nodes],
#         model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i,k] - model[:δ][j,k]) * prerequisites_data.a_l[i,j,k] ≤ (1 - model[:z][branch_dictionary[Set([i, j])],k]) * M)

#     JuMP.@constraint(model, pl_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in switched_Transmission_nodes],
#         model[:pij][i, j, k] - Sbase * (B[i, j]) * (model[:δ][i,k] - model[:δ][j,k]) * prerequisites_data.a_l[i,j,k] ≥ -(1 - model[:z][branch_dictionary[Set([i, j])],k]) * M)

#     if max_ts !== Inf
#         JuMP.@constraint(model, maximum_allowed_switching_actions[k in prerequisites_data.k], sum(1 - model[:z][branch_dictionary[l],k] for l in switched_Transmission_nodes) ≤ max_ts)
#     end

#     aux_bus_set = prerequisites_data.aux_bus_set
#     Reconf_dict = prerequisites_data.reconf_dict
#     Nodes_set = prerequisites_data.nodes_set
#     Reconf_nodes = prerequisites_data.reconf_nodes
#     Coupler_nodes = prerequisites_data.coupler_nodes
#     Coupler_dict = prerequisites_data.coupler_dict
#     Sbase = prerequisites_data.base_MVA
#     default_off_reconf = prerequisites_data.default_off_reconf
#     Coupler_set = prerequisites_data.coupler_set

#     # 2.5 Switching constraint to avoid connecting an element to two busbars at the same time
#     JuMP.@constraint(model, no_circular_path_constraint[bus in aux_bus_set, k in prerequisites_data.k], sum(model[:z_l][l, k] for l in grid.Buses[bus].ConnectedLinesIDs if grid.Branches[l].BranchType == 1) == 1)
#     M_δ = 2 * π

#     # 2.6.1 Phase angle constraints accross all reconfiguration lines to be the same if the switch is 1
#     JuMP.@constraint(model, phase_angle_equivalence_l_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Reconf_nodes],
#         model[:δ][j, k] - M_δ * (1 - model[:z_l][Reconf_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

#     JuMP.@constraint(model, phase_angle_equivalence_l_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Reconf_nodes],
#         model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - model[:z_l][Reconf_dict[Set([i, j])], k]))

#     # 2.6.2 Phase angle constraints accross all couplers to be the same if the switch is 1
#     JuMP.@constraint(model, phase_angle_equivalence_c_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Coupler_nodes],
#         model[:δ][j, k] - M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]) ≤ model[:δ][i, k])

#     JuMP.@constraint(model, phase_angle_equivalence_c_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Coupler_nodes],
#         model[:δ][i, k] ≤ model[:δ][j, k] + M_δ * (1 - model[:z_c][Coupler_dict[Set([i, j])], k]))

#     # 2.7.1 Reconfiguration line capacity
#     JuMP.@constraint(model, reconf_cap_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Reconf_nodes],
#         -model[:z_l][Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

#     JuMP.@constraint(model, reconf_cap_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Reconf_nodes],
#         model[:pij][i, j, k] ≤ model[:z_l][Reconf_dict[Set([i, j])], k] * Sbase * grid.Branches[Reconf_dict[Set([i, j])]].rating)

#     # 2.7.2 Reconfiguration line capacity
#     JuMP.@constraint(model, coupler_cap_1[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Coupler_nodes],
#         -model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating ≤ model[:pij][i, j, k])

#     JuMP.@constraint(model, coupler_cap_2[i in Nodes_set, j in Nodes_set, k in prerequisites_data.k; Set([i, j]) in Coupler_nodes],
#         model[:pij][i, j, k] ≤ model[:z_c][Coupler_dict[Set([i, j])], k] * Sbase * grid.Branches[Coupler_dict[Set([i, j])]].rating)

#     if max_reconf !== Inf
#         JuMP.@constraint(model, max_reconf_actions[k in prerequisites_data.k], sum(model[:z_l][r,k] for r in default_off_reconf) ≤ max_reconf)
#     end

#     if max_splitting !== Inf
#         JuMP.@constraint(model, max_splitting_actions[k in prerequisites_data.k], sum(model[:z_c][c,k] for c in Coupler_set) ≤ max_splitting)
#     end


#     # Load shedding limits based on bids
#     JuMP.@constraint(model, load_shedding_limits[d in load_set, k in prerequisites_data.k; grid.Loads[d].ActiveControl == true], 0 ≤ model[:p_ls][d,k] ≤ order_book.Load_bids[d]["qty"][t])
#     if no_load_shedding_pre_contingency
#         JuMP.@constraint(model, no_load_shedding_pre_contingency[d in load_set; grid.Loads[d].ActiveControl == true], model[:p_ls][d,1] == 0)
#     end
#     # GENERATOR CAPACITY
#     JuMP.@constraint(model, gen_active_power_limits_up[g in up_gen_set, k in prerequisites_data.k], 0 ≤ model[:p][g,k] ≤ order_book.Gen_bids[gen_id_dict[g]]["qty"][t][1] * prerequisites_data.a_g[g,k])
#     JuMP.@constraint(model, gen_active_power_limits_down[g in down_gen_set, k in prerequisites_data.k], order_book.Gen_bids[gen_id_dict[g]]["qty"][t][2] * prerequisites_data.a_g[g,k] ≤ model[:p][g,k] ≤ 0 )
#     # JuMP.@constraint(model, total_balance[k in prerequisites_data.k], sum(model[:p][g,k] for g in Gen_set) + sum(schedule["gen"][g][t] for g in Gen_set) == sum(grid.Loads[l].Pd for l in prereq.load_set))

#     Nodes_set = prerequisites_data.nodes_set
#     JuMP.@constraint(model, angle_limits[i in Nodes_set, k in prerequisites_data.k], grid.Buses[i].δ_min ≤ model[:δ][i,k] ≤ grid.Buses[i].δ_max)
#     if !isnothing(prerequisites_data.reference_node)
#         JuMP.@constraint(model, reference_node[k in prerequisites_data.k], model[:δ][prerequisites_data.reference_node,k] == 0)
#     end

#     non_base_contingencies = collect(setdiff(Set(prerequisites_data.k),Set([1])))

#     JuMP.@constraint(model,no_redispatch[g in prerequisites_data.gen_set, k in non_base_contingencies], model[:p][g,1] == model[:p][g,k])

#     # Corrective transmission switching constraints
#     if simulation_settings.transmission_switching == [:post]
#         JuMP.@constraint(model, static_grid_pre_contingency[l in prerequisites_data.switched_transmission_set], model[:z][l,1] == 1)
#     elseif simulation_settings.transmission_switching == [:pre]
#         JuMP.@constraint(model, static_grid_post_contingency[l in prerequisites_data.switched_transmission_set,k in non_base_contingencies], model[:z][l,k] == model[:z][l,1])
#     end

#     # Corrective substation switching constraints
#     if simulation_settings.substation_switching["reconf"] == [:pre]
#         JuMP.@constraint(model, static_substation_post_contingency[l in prerequisites_data.reconf_set, k in non_base_contingencies], model[:z_l][l,k] == model[:z_l][l,1])
#     end

#     if simulation_settings.substation_switching["splitting"] == [:pre]
#         JuMP.@constraint(model, static_coupler_post_contingency[l in prerequisites_data.coupler_set, k in non_base_contingencies], model[:z_c][l,k] == model[:z_c][l,1])
#     elseif simulation_settings.substation_switching["splitting"] == [:post]
#         JuMP.@constraint(model, static_coupler_pre_contingency[l in prerequisites_data.coupler_set], model[:z_c][l,1] == 1)
#     elseif simulation_settings.substation_switching["splitting"] == []
#         JuMP.@constraint(model, static_coupler[l in prerequisites_data.coupler_set, k in prerequisites_data.k], model[:z_c][l,k] == 1)
#     end

#     JuMP.@objective(model,Min, sum(model[:p][g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in up_gen_set) 
#         + sum(-model[:p][g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in down_gen_set) 
#         + sum(model[:p_ls][d,k]*order_book.Load_bids[d]["price"][t] for d in load_set , k in prerequisites_data.k if grid.Loads[d].ActiveControl == true))

#     JuMP.optimize!(model)

#     #===================================== Extracting Results =====================================#
#     if "Model" in keys(order_book.Grid_Results)
#         push!(order_book.Grid_Options["Model"], t => model)
#     else
#         order_book.Grid_Options["Model"] = Dict(t => model)
#     end
    
#     ΔPg = JuMP.value.(model[:p])
#     ΔPd = JuMP.value.(model[:p_ls])
    
#     ΣΔPg = sum(ΔPg[g,1] for g in prerequisites_data.gen_set)
#     if abs(ΣΔPg) ≥ 1e-4
#         error("Total amount of redispatch must add up to zero! but added up to -> $ΣΔPg")
#     end

#     ΣΔPd = Dict()
#     is_secure = 1
#     for k in prerequisites_data.k
#         ΣΔPd_k = sum(ΔPd[d,k] for d in load_set)
#         if abs(ΣΔPd_k) > 1e-4
#             is_secure = 0
#             push!(ΣΔPd, k => ΣΔPd_k)
#         end
#     end
    
#     SecurityStatus = Dict("Secure" => is_secure, "CriticalContingencies" => ΣΔPd)
#     if "SecurityStatus" in keys(order_book.Grid_Results)
#         push!(order_book.Grid_Results["SecurityStatus"], t => SecurityStatus)
#     else
#         order_book.Grid_Results["SecurityStatus"] = Dict(t => SecurityStatus)
#     end

#     ΣPd = sum(schedule["load"][d][t] + ΔPd[d,1] for d in load_set)

#     opex = sum(ΔPg[g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in up_gen_set) 
#         + sum(-ΔPg[g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in down_gen_set)
#         + sum(ΔPd[d,1]*order_book.Load_bids[d]["price"][t] for d in load_set if grid.Loads[d].ActiveControl == true)
    
#     unshed_opex = sum(ΔPg[g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in up_gen_set) 
#         + sum(-ΔPg[g,1]*order_book.Gen_bids[gen_id_dict[g]]["price"][t] for g in down_gen_set)

#     MCP = unshed_opex/ΣPd

#     lmp_i = Dict()
#     if has_duals(model)
#         for i in Nodes_set
#             push!(lmp_i, i => sum(dual(model[:Pnodal][i, k]) for k in prerequisites_data.k))
#         end
#     else
#         for i in Nodes_set
#             push!(lmp_i, i => MCP)
#         end
#     end

#     ########## Base Schedule Adjustment ##########
#     for duplet_id in keys(grid.Generator_Duplets)
#         gen_id = grid.Generator_Duplets[duplet_id][1]
#         gen_id_prime = grid.Generator_Duplets[duplet_id][2]
#         push!(order_book.Schedule["gen"][gen_id], t => ΔPg[gen_id,1] + ΔPg[gen_id_prime,1])
#     end

#     for load_id in order_book.Load_ids
#         push!(order_book.Schedule["load"][load_id], t => ΔPd[load_id,1])
#     end

#     return MCP, opex, lmp_i
# end