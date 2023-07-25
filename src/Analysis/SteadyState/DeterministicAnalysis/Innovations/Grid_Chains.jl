function calculate_line_profits(grid::PowerGrid, prerequisites; t=1, k=1, real_duals=true)
    update_grid_tables_DOPF!(grid; t=t, k=k)
    data = DataFrame(line_id=Int64[],FR=Int64[],TO=Int64[], λ_FR=Float64[], λ_TO=Float64[],λ_L=Float64[],
        Δλ=Float64[],PF=Float64[],Profit=Float64[],Utilization=Float64[])

    for row in eachrow(grid.LineLoading)
        
        branch_id = row[:BranchID]
        if grid.Branches[branch_id].BranchType == 0
            PF = row[:PL_fr_to]
            Utilization = row[:Utilization]
            if PF ≥ 0
                fr = row[:FromBus]
                to = row[:ToBus]
            else
                fr = row[:ToBus]
                to = row[:FromBus]
            end

            PF = abs(PF)
            
            if real_duals
                λ_FR = grid.Bus_Duals[t][fr]
                λ_TO = grid.Bus_Duals[t][to]
                λ_L  = grid.Line_Duals[t][branch_id]
            else
                λ_FR = prerequisites.Order_Book.λ_i_t[t][fr]
                λ_TO = prerequisites.Order_Book.λ_i_t[t][to]
                λ_L  = 0
            end
            Δλ   = λ_TO - λ_FR
            Profit = Δλ*PF
            push!(data, [branch_id, fr, to, λ_FR, λ_TO, λ_L, Δλ, PF, Profit, Utilization])
        end
    end
    return data
end

function calculate_line_profits_from_PF(grid::PowerGrid, power_flow_dict::Dict, profits_pre_switching::DataFrame)
    
    data = DataFrame(line_id=Int64[],FR=Int64[],TO=Int64[], λ_FR=Float64[], λ_TO=Float64[],
    Δλ=Float64[],PF=Float64[],Profit=Float64[],Utilization=Float64[],
    ΔPF=Float64[], ΔProfit=Float64[], ΔUtiliziation=Float64[], IsFlowFlipped=Float64[])

    for row in eachrow(profits_pre_switching)
        branch_id = row[:line_id]
        flow = power_flow_dict[branch_id][1]
        utilization = abs(power_flow_dict[branch_id][2])

        fr_ = splitting_grid.Branches[branch_id].Fr_bus_ID
        to_ = splitting_grid.Branches[branch_id].To_bus_ID
        change_dir_flag = 0
        if flow ≥ 0
            fr = fr_
            to = to_
            if row[:FR] == fr
                λ_FR = row[:λ_FR]
                λ_TO = row[:λ_TO]
            else
                λ_FR = row[:λ_TO]
                λ_TO = row[:λ_FR]
                change_dir_flag = 1
            end
        else
            fr = to_
            to = fr_
            if row[:FR] == fr
                λ_FR = row[:λ_FR]
                λ_TO = row[:λ_TO]
            else
                λ_FR = row[:λ_TO]
                λ_TO = row[:λ_FR]
                change_dir_flag = 1
            end
        end

        PF = abs(flow)
        Δλ   = λ_TO - λ_FR
        Profit = Δλ*PF

        if change_dir_flag == 1
            ΔPF = PF + row[:PF]
        else
            ΔPF = PF - row[:PF]
        end

        ΔProfit = Profit - row[:Profit]
        ΔUtiliziation = utilization - row[:Utilization]

        push!(data, [branch_id, fr, to, λ_FR, λ_TO, Δλ, PF, Profit, utilization, ΔPF, ΔProfit, ΔUtiliziation, change_dir_flag])
    end
    return data
end

function analyze_grid_chains(grid::PowerGrid, Profit_DF::DataFrame)
    """
    In this context, a `cycle` is a collection of nodes that form a closed loop, while a `chain` is the collection of corresponding branches
    in a `cycle`.
    """
    g = Graph(grid.N_bus)
    visited = []
    edge_labels = []

    for row in eachrow(Profit_DF)
        fr = row[:FR]
        to = row[:TO]
        pf = round(row[:PF],digits=2)
        U = round(row[:Utilization],digits=2)
        profit = round(row[:Profit],digits=2)
        if (fr,to) ∉ visited
            add_edge!(g, fr, to)
            my_tuple = (pf,profit,U)
            push!(edge_labels, "$my_tuple")
            push!(visited, (fr,to))
        end
    end

    cycles = cycle_basis(g)

    line_dict = Dict()
    for line_id in keys(grid.Branches)
        fr = grid.Branches[line_id].Fr_bus_ID
        to = grid.Branches[line_id].To_bus_ID
        key_set = Set([fr, to])
        if key_set ∉ keys(line_dict)
            push!(line_dict, key_set => [line_id])
        else
            push!(line_dict[key_set], line_id)
        end
    end

    profit_in_cycle = []
    Δλ_in_cycle = []
    PF_in_cycle = []
    chains = Dict()
    chain_id = 0
    for cycle in cycles
        chain_id += 1
        lines_in_cycle = []
        for i in range(1,length(cycle))
            if i == length(cycle)
                append!(lines_in_cycle, line_dict[Set([cycle[i], cycle[1]])])
            else
                append!(lines_in_cycle, line_dict[Set([cycle[i], cycle[i+1]])])
            end
        end

        
        Σprofit = 0
        ΣΔλ = 0
        ΣPF = 0
        for line_id in lines_in_cycle
            df_filtered = filter(row -> row.line_id == line_id, Profit_DF)
            Σprofit += df_filtered[1,:Profit]
            ΣΔλ += df_filtered[1,:Δλ]
            ΣPF += df_filtered[1,:PF]
        end

        push!(profit_in_cycle, Σprofit)
        push!(Δλ_in_cycle, ΣΔλ)
        push!(PF_in_cycle, ΣPF)
        push!(chains, chain_id => (lines_in_cycle, cycle, Σprofit, ΣΔλ, ΣPF))
    end

    return chains, cycles, profit_in_cycle, Δλ_in_cycle, PF_in_cycle
end

function get_sub_graphs(grid::PowerGrid, chains::Dict)
    """
    A `sub-graph` about node `n` is defined as all chains that intersect at least at node `n`.
    Output: Dictionary where := key (node id) => value (vector of chain ids intersect at "node id") 
    """
    sub_graphs = Dict()
    for node_id in keys(grid.Buses)
        push!(sub_graphs, node_id => [])
    end

    for node_id in keys(grid.Buses)
        for chain_id in keys(chains)
            cycle = chains[chain_id][2]
            if node_id ∈ cycle
                append!(sub_graphs[node_id], chain_id)
            end
        end
    end

    return sub_graphs
end