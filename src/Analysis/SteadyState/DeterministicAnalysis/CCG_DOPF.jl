# include("DOPF_Constraints.jl")

function build_DOPF_MP!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites)
    return build_full_DOPF_model!(grid,SimulationSettings,prerequisites_data)
end

function build_DOPF_SP!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    solved_MP_model::Model, t::Int, k::Int)
    prerequisites_data_instance = deepcopy(prerequisites_data)
    prerequisites_data_instance.time_horizon = [t]
    prerequisites_data_instance.k_t[t] = [k]

    # # Fixing schedules and emptying commitable/non-commitable gen lists:
    # p_gen_ac = JuMP.value.(solved_MP_model[:p_gen_ac][:,1,t])
    # fixed_schedules = Dict()
    # for g in prerequisites_data.ac_gen_ids
    #     push!(fixed_schedules, )
    # # Fixing relevant topology:
    # if SimulationSettings.

    # Building model:

end

function solve_decomposed_model!(model::Model)
    optimize!(model)
    return model
end

function solve_DOPF_CCG!(grid::PowerGrid, SimulationSettings::DOPF_SimulationSettings, prerequisites_data::DOPF_Prerequisites,
    order_book::OrderBook ; update_grid=false, ϵ=0.1)

    δ = Inf
    
    iter_count = 0
    T = SimulationSettings.time_horizon
    K_all = prerequisites_data.k

    t_master = []
    t_sub = []

    δ_it = []
    LS_it = [] # to be calculated later

    solved_MP_model = []
    #Main loop
    while δ ≥ ϵ
        iter_count = iter_count + 1

        println("Debug msg @ iteration " * string(iter_count))
        # Solve master-problem

        MP_model = build_DOPF_MP!(grid, SimulationSettings, prerequisites_data)
        t_master_now = @elapsed solved_MP_model = solve_decomposed_model!(MP_model)
        push!(t_master, t_master_now)

        LB_t = [calculate_opex_t(solved_MP_model,grid, prerequisites_data, t; include_load_shedding=true) for t in T]

        UB_tk = zeros(length(T), length(K_all))

        t_sub_now = @elapsed Threads.@threads for t in T
            Threads.@threads for k in collect(setdiff(Set(K_all),Set(prerequisites_data.k_t[t])))
                SP_model = build_DOPF_SP!(grid, SimulationSettings, prerequisites_data)
                solved_SP_model = solve_decomposed_model!(SP_model)
                UB_tk[t,k] = JuMP.objective_value(solved_SP_model)
            end
        end
        push!(t_sub, t_sub_now)
        UB_t = zeros(length(T),1)
        
        for t in T
            # 1 should be a data structure ::CCG_Settings which can include nk[t] (number of k's to include) or a DRL agent that chooses the k's
            append!(prerequisites_data.k_t[t], find_next_k(UB_tk[t,:], 1))
            UB_t[t] = maximum(UB_tk[t,:])
        end

        δ = maximum(abs.(UB_t - LB_t) ./ LB_t)
        push!(δ_it, δ)
        
        println("|Iteration Count|δ|")
        println(iter_count, "  ", δ)
    end

    return solved_MP_model
end


function find_next_k(v::Vector, nk::Int)
    return sortperm(v, rev=true)[1:nk]
end
