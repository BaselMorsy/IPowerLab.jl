function calculate_powerflow(grid::PowerGrid; model=:DC, method=:NR, t=1, k=1, ref_node=1)
    δ = zeros(grid.N_bus)
    V = ones(grid.N_bus)
    Pl = Dict()
    Ql = Dict()
    if model == :DC

        B = B_Bus_Grid(grid)
        Pd = zeros(grid.N_bus, 1)
        for d in sort(collect(keys(grid.Loads)))
            load_bus = grid.Loads[d].LoadBus_ID
            Pd[load_bus] = grid.Loads[d].Pd_t[t]
        end
        
        Pg = zeros(grid.N_bus, 1)
        for g in sort(collect(keys(grid.Generators)))
            Pg[grid.Generators[g].GenBus_ID] = grid.Generators[g].Pg_tk[t][k]
        end

        P_inj = Pg - Pd
        θ = -inv(B[1:end .!= ref_node, 1:end .!= ref_node])*(P_inj[1:end .!= ref_node]./grid.S_base)
        insert!(θ, ref_node, 0)

        for l in keys(grid.Branches)
            fr = grid.Branches[l].Fr_bus_ID
            to = grid.Branches[l].To_bus_ID
            P_flow = grid.S_base*(1/grid.Branches[l].x)*(θ[fr] - θ[to])
            utilization = 100*(P_flow/(grid.S_base*grid.Branches[l].rating))
            push!(Pl, l => (P_flow, utilization))
        end

        return θ, Pl
    elseif model == :AC

    elseif model == :FD

    else
        error("Invalid power flow model ($model) !")
    end
    
end

