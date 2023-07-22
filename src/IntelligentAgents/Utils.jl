abstract type Agent end
# Define a single experience
struct Experience
    # State, action, reward, next state, done flag
    state::Array{Float64, 1}
    action::Array{Float64, 1}
    reward::Float64
    next_state::Array{Float64, 1}
    done
end

# Define the replay buffer
mutable struct ReplayBuffer
    # Maximum number of experiences to store
    capacity::Int
    # Stored experiences
    experiences::Array{Experience, 1}
    # Current index to store the next experience
    index::Int
    # Number of stored experiences
    size::Int
end

# Function to store an experience in the replay buffer
function store!(buffer::ReplayBuffer, experience::Experience)

    if buffer.size < buffer.capacity
        push!(buffer.experiences,experience)
        buffer.index += 1
        buffer.size += 1
    else
        # Store the experience at the current index
        buffer.experiences[buffer.index] = experience
        # Update the current index and size
        buffer.index = mod1(buffer.index + 1, buffer.capacity)
        buffer.size = min(buffer.size + 1, buffer.capacity)
    end
end

# Function to store an experience in the replay buffer
function store_experience!(agent::Agent, state::Array{Float64, 1}, action::Array{Float64, 1}, reward::Float64, next_state::Array{Float64, 1}, done)
    # Create an experience from the input data
    experience = Experience(state, action, reward, next_state, done)
    # Store the experience in the replay buffer
    store!(agent.replay_buffer, experience)
end

# Function to initialize the replay buffer
function ReplayBuffer(capacity::Int)
    # Create an empty replay buffer
    return ReplayBuffer(capacity, Experience[], 1, 0)
end

# Function to sample a batch of experiences from the replay buffer
function sample(buffer::ReplayBuffer, batch_size::Int)
    # Sample a batch of experiences from the replay buffer
    indices = rand(1:buffer.size, batch_size)
    return buffer.experiences[indices]
end

