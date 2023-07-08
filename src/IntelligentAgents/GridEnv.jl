
abstract type Environment end

using Parameters
@with_kw mutable struct OPF_ENV <: Environment
    grid::PowerGrid
    prerequisites::OPF_Prerequisites
    simulation_settings::OPF_SimulationSettings
    n_states::Int64
    n_actions::Int64
end

function OPF_ENV_init!(grid::PowerGrid, agent_type::Symbol, prerequisites::OPF_Prerequisites, simulation_settings::OPF_SimulationSettings)

    if agent_type == :DDPG
        state_size = 2 * grid.N_bus
        action_size = Int(length(prerequisites.Coupler_set) + length(prerequisites.reconf_set) / 2)
        actor_learning_rate = 1e-3
        critic_learning_rate = 1e-3
        discount_factor = 0.99
        noise_theta = 0.1
        noise_sigma = 0.2
        replay_buffer_capacity = 10000
        agent = DDPGAgent(state_size, action_size, actor_learning_rate, critic_learning_rate,
            discount_factor, noise_theta, noise_sigma, replay_buffer_capacity)
        my_env = OPF_ENV(grid, prerequisites, simulation_settings, state_size, action_size)

        return my_env, agent
    elseif agent_type == :DQN
    end
end

function get_state(env::Environment)
    state = zeros(env.grid.N_bus * 2)
    for bus_id in sort(collect(keys(env.grid.Buses)))
        gen_ids = env.grid.Buses[bus_id].ConnectedGensIDs
        load_ids = env.grid.Buses[bus_id].ConnectedLoadsIDs
        if load_ids == []
            state[bus_id] = 0
        else
            state[bus_id] = sum(env.grid.Loads[l].Pd for l in load_ids)
        end
        if gen_ids == []
            state[bus_id+env.grid.N_bus] = 0
        else
            state[bus_id+env.grid.N_bus] = sum(env.grid.Generators[g].Pg for g in gen_ids)
        end
    end
    return state
end

function get_state(env::Environment,pg_dict)
    state = zeros(env.grid.N_bus * 2)
    for bus_id in sort(collect(keys(env.grid.Buses)))
        gen_ids = env.grid.Buses[bus_id].ConnectedGensIDs
        load_ids = env.grid.Buses[bus_id].ConnectedLoadsIDs
        if load_ids == []
            state[bus_id] = 0
        else
            state[bus_id] = sum(env.grid.Loads[l].Pd for l in load_ids)
        end
        if gen_ids == []
            state[bus_id+env.grid.N_bus] = 0
        else
            state[bus_id+env.grid.N_bus] = sum(get(pg_dict,g,0) for g in gen_ids)
        end
    end
    return state
end

function extract_switching_from_action(action, prerequisites_data, grid)
    z_l = Dict()
    z_c = Dict()

    twins = Dict()

    i = 1
    for substation_id in sort(collect(keys(grid.Substations)))
        push!(z_c, grid.Substations[substation_id].Reconf_CouplerLines_IDs[1] => round(action[i]))
        i += 1
        merge!(twins, grid.Substations[substation_id].twins)
    end

    for line_id in sort(prerequisites_data.default_off_reconf)
        push!(z_l, line_id => round(action[i]))
        push!(z_l, twins[line_id] => 1 - round(action[i]))
        i += 1
    end

    return z_l, z_c

end

function solve_static_topology_OPF!(grid, action, simulation_settings, prerequisites_data)

    z_l, z_c = extract_switching_from_action(action, prerequisites_data, grid)
 
    # Model initialization
    model = Model(simulation_settings.MILP_solver)

    # Variables initialization
    single_period_OPF_variable_initialization!(model, simulation_settings, prerequisites_data)

    # AC grid constraints:
    single_period_angle_limits_ac_grid!(prerequisites_data, grid, model)
    single_period_voltage_limits_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_generator_limits_ac_grid!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_nodal_balance_ac_node!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    single_period_transmission_capacity_limits_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)
    single_period_powerflow_ac_branch!(prerequisites_data, grid, model, simulation_settings.ac_grid_model)

    single_period_reconf_split_SG_constraints_ac_grid!(prerequisites_data, grid, model, z_l, z_c)

    # DC grid constraints:
    if length(keys(grid.DCBuses)) !== 0
        single_period_voltage_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_powerflow_dc_branch!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_transmission_capacity_limits_dc_grid!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
        single_period_nodal_balance_dc_node!(prerequisites_data, grid, model, simulation_settings.dc_grid_model)
    end

    # Converter constraints:
    if length(keys(grid.Converters)) !== 0 || length(keys(grid.DCLinks)) !== 0
        single_period_converter_constraints!(prerequisites_data, grid, model, simulation_settings.converter_model)
    end

    single_period_objective!(prerequisites_data, grid, model, simulation_settings.transmission_switching, simulation_settings.substation_switching)


    optimize!(model)
    if JuMP.has_values(model)
        _ = OPF_post_processing!(grid, model, simulation_settings, prerequisites_data)
        update_grid_tables!(grid)
    end

    return model
end

function step!(env::Environment, action::Array{Float64,1})

    model = solve_static_topology_OPF!(env.grid, action, env.simulation_settings, env.prerequisites)
    reward = 0
    done = 0
    if !JuMP.has_values(model)
        done = 1
        reward = -sum(env.grid.Generators[g].C1 * env.grid.Generators[g].Pg_max for g in env.prerequisites.Gen_set)
    else
        done = 0
        reward = -JuMP.objective_value(model)
    end

    state = get_state(env)
    return state, reward, done
end

function run!(agent,env,episodes,mode=:train;Pg_dict=nothing)
    for e in collect(1:episodes)

        done = 0
        reset_noise!(agent)
        k = 0
        while done == 0
            k += 1
            if isnothing(Pg_dict)
                state = get_state(env)
            elseif isa(Pg_dict,Dict)
                state = get_state(env, Pg_dict)
            end

            action = select_action(agent,state,mode=mode)
            reward_prev_index = agent.replay_buffer.index - 1
            next_state, reward, done = step!(env, action)

            if abs(reward - agent.replay_buffer.experiences[reward_prev_index].reward)/reward ≤ 0.05
                done = 1
            end

            if k ≥ 5
                done = 1
            end

            store_experience!(agent, state, action, reward, next_state, done)
            if mode == :train
                train!(agent, 25)
            end

            if isa(Pg_dict, Dict) && done == 1
                return extract_switching_from_action(action, env.prerequisites, env.grid)
            end
        end
    end
end

function run_single!(agent, env, mode; Pg_dict=nothing)
    done = 0
    # reset_noise!(agent)
    if isnothing(Pg_dict)
        state = get_state(env)
    elseif isa(Pg_dict, Dict)
        state = get_state(env, Pg_dict)
    end

    action = select_action(agent, state, mode=mode)
    reward_prev_index = agent.replay_buffer.index - 1
    f = state - agent.replay_buffer.experiences[reward_prev_index].state
    indexes = f .!= 0
    different_state_vars = f[indexes]
    gen_prod = state[119:end]
    if sum(different_state_vars) / length(different_state_vars) ≤ 0.05 * sum(gen_prod)
        z_l, z_c = extract_switching_from_action(action, env.prerequisites, env.grid)
        return z_l, z_c, 1
    end
    next_state, reward, done = step!(env, action)

    if abs(reward - agent.replay_buffer.experiences[reward_prev_index].reward) / abs(reward) ≤ 0.05
        done = 1
    end

    store_experience!(agent, state, action, reward, next_state, done)
    if mode == :train
        train!(agent, 25)
    end
    z_l, z_c = extract_switching_from_action(action, env.prerequisites, env.grid)
    return z_l, z_c, 0
end

function solve_assisted_OPF!(env, agent)
    state = get_state(env)
    z_l, z_c = extract_switching_from_action(select_action(agent, state, mode=:test), env.prerequisites, env.grid)
    model = build_OPF_model!(env.grid, env.simulation_settings, env.prerequisites, optimize=false, clean=false)
    for coupler in env.prerequisites.Coupler_set
        set_start_value(model[:z_c][coupler], z_c[coupler])
    end
    for reconf in env.prerequisites.reconf_set
        set_start_value(model[:z_l][reconf], z_l[reconf])
    end
    optimize!(model)

    return model
end