using BSON
mutable struct OrnsteinUhlenbeckProcess
    # Number of dimensions
    dimensions::Int
    # Mean reversion rate
    theta::Float64
    # Standard deviation
    sigma::Float64
    # Current state
    state::Array{Float64,1}
end

# Define the DDPG agent
mutable struct DDPGAgent <: Agent
    # Network for the actor
    actor_network::Chain
    # Network for the critic
    critic_network::Chain
    # Target networks for the actor and critic
    actor_target_network::Chain
    critic_target_network::Chain
    # Learning rate for the actor and critic networks
    actor_learning_rate::Float64
    critic_learning_rate::Float64
    # Discount factor
    discount_factor::Float64
    # Noise process for exploration
    noise_process::OrnsteinUhlenbeckProcess
    # Replay buffer for storing past experiences
    replay_buffer::ReplayBuffer
    actor_opt
    critic_opt
end

# Function to create the actor network
function create_actor_network(state_size::Int, action_size::Int)
    # Define the network layers
    layers = [
        Dense(state_size, 32, relu),
        Dense(32, 32, relu),
        Dense(32, action_size, sigmoid)
    ]
    # Create and return the network
    return Chain(layers...)
end

# Function to create the critic network
function create_critic_network(state_size::Int, action_size::Int)
    # Define the network layers
    layers = [
        Dense(state_size + action_size, 32, relu),
        Dense(32, 32, relu),
        Dense(32, 1)
    ]
    # Create and return the network
    return Chain(layers...)
end

# Initialize the Ornstein-Uhlenbeck process
function OrnsteinUhlenbeckProcess(dimensions::Int, theta::Float64, sigma::Float64)
    # Initialize the state to 0
    state = zeros(dimensions)
    # Return the Ornstein-Uhlenbeck process
    return OrnsteinUhlenbeckProcess(dimensions, theta, sigma, state)
end

# Generate a sample from the process
function (process::OrnsteinUhlenbeckProcess)()
    # Update the state
    process.state = process.state .+ process.theta .* (0 .- process.state) .+ process.sigma .* randn(process.dimensions)
    # Return the state
    return process.state
end

# Reset the process
function reset!(process::OrnsteinUhlenbeckProcess)
    # Reset the state to 0
    process.state .= 0
end

# Initialize the DDPG agent
function DDPGAgent(state_size::Int, action_size::Int, actor_learning_rate::Float64, critic_learning_rate::Float64, discount_factor::Float64, noise_theta::Float64, noise_sigma::Float64, replay_buffer_capacity::Int)
    # Create the actor and critic networks
    actor_network = create_actor_network(state_size, action_size)
    critic_network = create_critic_network(state_size, action_size)
    # Create the target networks
    actor_target_network = deepcopy(actor_network)
    critic_target_network = deepcopy(critic_network)
    # Create the noise process
    noise_process = OrnsteinUhlenbeckProcess(action_size, noise_theta, noise_sigma)
    # Create the replay buffer
    replay_buffer = ReplayBuffer(replay_buffer_capacity)
    actor_opt = Adam(0.001, (0.9, 0.8))
    critic_opt = Adam(0.001, (0.9, 0.8))
    # Return the DDPG agent
    return DDPGAgent(actor_network, critic_network, actor_target_network, critic_target_network, actor_learning_rate, critic_learning_rate, discount_factor, noise_process, replay_buffer, actor_opt, critic_opt)
end

# Function to select an action using the actor network
function select_action(agent::DDPGAgent, state::Array{Float64,1}; mode=:train)
    # Predict the action using the actor network
    action = agent.actor_network(state)
    if mode == :train
        # Add noise for exploration
        action += agent.noise_process()
    end
    # Clamp the action within the valid range
    action = round.(clamp.(action, 0, 1.0))
    # Return the action
    return action
end

# Function to train the DDPG agent
function train!(agent::DDPGAgent, batch_size::Int)
    # Check if there are enough experiences in the replay buffer
    if agent.replay_buffer.size < batch_size
        return
    end
    # Sample a batch of experiences from the replay buffer
    experiences = sample(agent.replay_buffer, batch_size)
    # Convert the experiences to arrays
    states = [e.state for e in experiences]
    actions = [e.action for e in experiences]
    rewards = [e.reward for e in experiences]
    next_states = [e.next_state for e in experiences]
    dones = [e.done for e in experiences]
    # Predict the next actions using the target actor network
    next_actions = agent.actor_target_network.(next_states)
    # Predict the Q-values using the target critic network
    q_values = []
    q_values_actual = []
    for i in 1:length(next_states)
        q_input = vcat(next_states[i], next_actions[i])
        push!(q_values, agent.critic_target_network(q_input)[1])

        q_input_act = vcat(states[i], actions[i])
        push!(q_values_actual, agent.critic_network(q_input_act)[1])
    end

    # Calculate the expected Q-values
    expected_q_values = rewards .+ agent.discount_factor .* (1 .- dones) .* q_values

    # Calculate the critic loss
    # critic_loss = mean((expected_q_values .- q_values_actual) .^ 2)
    critic_loss_func(x, y) = mean((x .- y) .^ 2)

    # critic params
    critic_params = Flux.params(agent.critic_network)

    #critic data
    critic_data = [(expected_q_values[i], q_values_actual[i]) for i in 1:length(expected_q_values)]
    # update critic network
    Flux.@epochs 10 Flux.train!(critic_loss_func, critic_params, critic_data, agent.critic_opt)

    # Predict the actions using the actor network
    predicted_actions = agent.actor_network.(states)
    q_pred = []
    for i in 1:length(states)
        input = vcat(states[i], predicted_actions[i])
        push!(q_pred, agent.critic_network(input))
    end
    # Calculate the actor loss
    actor_loss_func(x) = -mean(x)
    actor_params = Flux.params(agent.actor_network)
    # update actor network
    Flux.@epochs 10 Flux.train!(actor_loss_func, actor_params, q_pred, agent.actor_opt)
    # Update the target networks
    update_target_networks!(agent.actor_target_network, agent.actor_network)
    update_target_networks!(agent.critic_target_network, agent.critic_network)
end

# Function to update a target network
function update_target_networks!(target_network::Chain, network::Chain)
    # Update the target network weights
    for (target_param, param) in zip(Flux.params(target_network), Flux.params(network))
        target_param .= 0.001 .* param .+ 0.999 .* target_param
    end
end

# Function to reset the noise process
function reset_noise!(agent::DDPGAgent)
    reset!(agent.noise_process)
end

function save_agent(agent::DDPGAgent, path)
    save(path*".jld2","ddpg",agent)
end

function load_agent(agent_path)
    return load(agent_path)["ddpg"]
end