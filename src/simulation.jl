include("model.jl") 

using JLD

mutable struct SimulationData
    step::Int32
    clusters::Int32
    size::Vector{Int32}
    gyr::Vector{Float64}
end

mutable struct EnsembleData
    mrna_ang::Float64
    n_lip::Int32
    sim_data::Array{SimulationData}
end

# Initializer for empty data types
SimulationData(step) = SimulationData(step, 0, Vector{Int32}(), Vector{Float64}())
EnsembleData(sim_data) = EnsembleData(0.0, 0, sim_data)


# Data collection functions for simulation
function numMembers(agent)
    # Determines the number of members in a cluster
    # A free liposome is considered to have a size of 1
    if agent.type == :Cluster
        return length(agent.mems)
    elseif agent.type == :Member
        return 0
    elseif agent.type == :Free
        return 1
    end
end
                
function gyr_radius(agent)
    # Returns the radius of gyration for clusters and 0 otherwise
    if agent.type == :Cluster
        return agent.info[1]
    else
        return 0.0
    end
end

agentData = [numMembers, gyr_radius] # array of functions to be passed into Agents.jl package

function _run!(numSteps; skip = 1, seed = 0, mrna_ang = pi/2, filename = "")
    model = initialize_model(n_liposomes = 125, seed = seed)
    
    data = Array{SimulationData}([SimulationData(step) for step in 1:skip:numSteps])
    
    for i in 1:Int(floor(numSteps/skip))
        step!(model, dummystep, system_step!, skip)
        for id in allids(model)
            if model[id].type == :Cluster
                data[i].clusters += 1
                push!(data[i].size, length(model[id].mems))
                push!(data[i].gyr, 100 * model[id].info[1])
            end
        end
    end

    if filename == ""
        @save joinpath(@__DIR__, "../data/single_run.jld") data
    else
        @save filename data
    end
    return data
end
    
function _run!(numSteps, data; skip = 1, seed = 0, mrna_ang = pi/2, n_liposomes = 125)
    model = initialize_model(n_liposomes = n_liposomes, seed = seed, mrna_ang = mrna_ang)
        
    for i in 1:Int(floor(numSteps/skip))
        step!(model, dummystep, system_step!, skip)
        for id in allids(model)
            if model[id].type == :Cluster
                data[i].clusters += 1
                push!(data[i].size, length(model[id].mems))
                push!(data[i].gyr, 100 * model[id].info[1])
            end
        end
    end
    return data
end
                
function _ensemble!(numTests::Integer, numSteps::Integer, mrna_angs::Real, n_liposomes::Integer; skip = 1, filename = "", iter = 1)
    n_liposomes = fill(n_liposomes, numTests)
    mrna_angs = fill(mrna_angs, numTests)
    _ensemble!(numTests, numSteps, mrna_angs, n_liposomes; skip=skip, filename=filename, iter=iter)
end
                
function _ensemble!(numTests::Integer, numSteps::Integer, mrna_angs::AbstractArray, n_liposomes::Integer; skip = 1, filename = "", iter = 1)
    length(mrna_angs) == numTests || error("Length of array of mRNA angles does not match number of runs")
    n_liposomes = fill(n_liposomes, numTests)
    _ensemble!(numTests, numSteps, mrna_angs, n_liposomes; skip=skip, filename=filename, iter=iter)
end
                
function _ensemble!(numTests::Integer, numSteps::Integer, mrna_angs::Real, n_liposomes::AbstractArray; skip = 1, filename = "", iter = 1)
    length(n_liposomes) == numTests || error("Number of liposomes array length does not match number of runs")
    mrna_angs = fill(mrna_angs, numTests)
    _ensemble!(numTests, numSteps, mrna_angs, n_liposomes; skip=skip, filename=filename, iter=iter)
end
                
function _ensemble!(numTests::Integer, numSteps::Integer, mrna_angs::AbstractArray, n_liposomes::AbstractArray; skip = 1, filename = "", iter = 1)
    # Input validation
    length(mrna_angs) == numTests || error("Length of array of mRNA angles does not match number of runs")
    length(n_liposomes) == numTests || error("Length of number of liposomes array does not match number of runs")

    seeds = rand(UInt8, numTests*iter) # Array of seeds for models
    empty_run = Array{SimulationData}([SimulationData(step) for step in 1:skip:numSteps]) # initialize data structure
    data = Array{EnsembleData}([EnsembleData(copy(empty_run)) for i in 1:numTests])

    for i in 1:numTests
        data[i].mrna_ang = mrna_angs[i]
        data[i].n_lip = n_liposomes[i]
        for j in 1:iter
            data[i].sim_data = _run!(numSteps, data[i].sim_data; skip = skip, seed = seeds[numTests * (j - 1) + i], mrna_ang = mrna_angs[i], n_liposomes = n_liposomes[i])
        end
    end
    
    if filename == ""
        @save joinpath(@__DIR__, "..\\data\\ensemble_run.jld") data
    else
        @save filename data
    end
    
    return data
end