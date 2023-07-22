struct DQNAgent <: Agent
    model::Chain
    target_model::Chain
    opt::ADAM
    γ::Float64  # discount factor
    replay_buffer::ReplayBuffer  # experience replay buffer
    epsilon::Float64  # probability of choosing a random action
    epsilon_min::Float64  # minimum value of epsilon
    epsilon_decay::Float64  # rate at which epsilon decays over time
    tau::Float64  # rate at which target model is updated
end

function DQNAgent(state_size::Int, action_size::Int;
                  γ=0.99, epsilon=1.0, epsilon_min=0.01, epsilon_decay=0.995,
                  tau=0.01,
                  buffer_size=10000, batch_size=64)
    layers = [
        Dense(state_size, 32, relu),
        Dense(32, 32, relu),
        Dense(32, action_size)
        ]
    # Create and return the network
    model = Chain(layers...)
    target_model = deepcopy(model)
    opt = Adam(0.001, (0.9, 0.8))
    exp_replay = ReplayBuffer(buffer_size)
    return DQNAgent(model, target_model, opt, γ, exp_replay, epsilon, epsilon_min, epsilon_decay, tau)
end

# function (agent::DQNAgent)(state)
#     if rand() < agent.epsilon
#         # choose a random action with probability epsilon
#         return rand(1:size(agent.model.layers[end].W, 2))
#     else
#         # choose the action with the highest predicted Q-value
#         return argmax(agent.model(state))
#     end
# end

# Function to select an action 
function (agent::DQNAgent)(state::Array{Float64, 1};train ::Bool = true)
    if train 
        if rand() < agent.epsilon
            # choose a random action with probability epsilon
            return rand(1:size(agent.model.layers[end].W, 2))
        else
            # choose the action with the highest predicted Q-value
            return argmax(agent.model(state))
        end
    else
        return argmax(agent.model(state))
    end
end

# function (agent::DQNAgent)(state, action, reward, next_state, done)
#     # store the experience in the replay buffer
#     agent.exp_replay(state, action, reward, next_state, done)
#     # sample a batch of experiences from the replay buffer
#     states, actions, rewards, next_states, dones = agent.exp_replay()
#     # compute the target Q-values for the sampled experiences
#     q_targets = agent.target_model(next_states)
#     for i in 1:size(q_targets, 2)
#         q_targets[:, i] = rewards[:, i] .+ agent.γ .* maximum(q_targets[:, i], 1 .- dones[:, i])
#     end
#     # compute the predicted Q-values for the sampled experiences
#     q_preds = agent.model(states)
#     for i in 1:size(q_preds, 2)
#         q_preds[:, i] = q_preds[:, i][actions[:, i]]
#     end
#     # compute the mean squared error between the predicted and target Q-values
#     loss = mse(q_preds, q_targets)
#     # update the model parameters to minimize the loss
#     Flux
# end