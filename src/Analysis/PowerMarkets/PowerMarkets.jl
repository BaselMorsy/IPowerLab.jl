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