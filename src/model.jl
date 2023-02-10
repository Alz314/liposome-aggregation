using Agents
using Random
using Rotations
using LinearAlgebra
using Distributions
using SpecialFunctions


const k_b = 1.380649e-23 # m^2 kg s^-2 K^-1 (Boltzmann constant)

mutable struct Liposome <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64} # velocity in x, y, and z directions
    type::Symbol # a liposome can either be a free liposome; a "cluster" object; or a member of a cluster
    
    mems::Vector{Int} # Multipurpose vector depending on whether this is a liposome or cluster
    # For clusters, mems is a vector of all member ids
    # For members, mems is a vector of the ids that it is specifically bound to
    info::Vector{Float64} 
    # For clusters, info contains the gyration radius, hydrodynamic radius, and rotational diffusion constant
    # For members, info contains the id of the cluster it belongs to and the number of members of that cluster
    oab::Vector{Float64}  # Orientation, mRNA point A, mRNA point B angles
end
    
# Constructors for each type of liposome
Free(id, pos, vel, oab) = Liposome(id, pos, vel, :Free, [], [], oab) # Free liposomes are not bound to any cluster
Cluster(id, pos, members) = Liposome(id, pos, (0.,0.), :Cluster, members, [0., 0., 0.], [0,0,0]) # Fake agent that acts as the CM of the cluster
Member(id, pos, cluster, oab) = Liposome(id, pos, (0.,0.), :Member, [], [cluster], oab) # Liposomes that are part of a cluster
    
get_noise(model) = (rand(model.rng, Normal(0,1)), rand(model.rng, Normal(0,1)))

function get_oab(model)
    oab = rand(3) * 2 * pi
    oab[3] = (oab[2] + model.mrna_ang) % (2*pi)
    return oab
end

function rand_vel(model)
    # Generates a uniformly distributed random unit vector
    θ = rand(model.rng) * 2π 
    vel = (cos(θ), sin(θ))
    return vel
end

function initialize_model(;
    n_liposomes = 100,
    speed = 1.0,
    extent = 100,
    spacing = 1,
    prob = 0.8,
    seed = 0,
    mass_liposome = 2.5e-18, 
    liposome_r = 100,
    viscosity = 0.01,
    T = 293,
    alpha = 2.0,
    mrna_ang = pi/2
)
    # If a seed is not set by the user, generate a random one
    rng = seed == 0 ? MersenneTwister(trunc(Int, rand()*10)) : MersenneTwister(seed)
    
    dt = 3 * π * viscosity * (liposome_r*1e-9)^3 / (k_b * T)
    rot_dif = k_b * T/ (8 * π * viscosity * (1e-7)^3 ) 
    #dt = 0.001
    
    # Setup of the model and space
    properties = (prob=prob, dt=dt, speed=speed, extent=extent, mass=mass_liposome, radius=liposome_r, T=T, eta=viscosity, alpha=alpha, rot_dif=rot_dif, mrna_ang=mrna_ang)
    space_extents = (extent, extent)
    scheduler = Schedulers.randomly
    space2d = ContinuousSpace(space_extents; periodic = true, spacing = spacing)
    model = ABM(Liposome, space2d; scheduler, properties, rng)
    
    # Adding in free liposomes using the Free() constructor
    for _ in 1:n_liposomes
        add_agent_pos!(Free(nextid(model), random_position(model), rand_vel(model), get_oab(model)), model)
    end
    
    return model
end

function free_liposome_step!(liposome, model)
    # Free liposome step first checks any interaction with other liposomes and then moves the liposome
    
    # Obtain the ids of neighbors that are close enough to interact with the liposome
    neighbor_ids = nearby_ids(liposome, model, 1)
    
    # Bind liposomes based on the type of the neighbor
    for id in neighbor_ids
        
        # Check to see if this interaction results in aggregation
        if model[id].type != :cluster
            if !successful_interaction(liposome.id, id, model)
                continue # if not, skip this neighbor
            end
        else
            continue
        end
        
        neighbor = model[id]
        
        # Binding the free liposome to an existing cluster
        if neighbor.type == :Member && liposome.type != :Member
            # Checking that the liposome isn't a member fixes a bug where a liposome can get counted twice
            clusterId = Int64(neighbor.info[1])
            cluster = model[clusterId]
            cluster.mems = [cluster.mems; liposome.id] # add the liposome to the cluster's list
            neighbor.mems = [neighbor.mems; liposome.id] # add the liposome to the neighbor's list
            liposome.type = :Member
            liposome.mems = [id]
            liposome.info = [clusterId, size(cluster.mems)[1]]
            adjustCM(cluster, model)
            break
        end
        
        # Binding two free liposomes
        if neighbor.type == :Free
            # Generate a new cluster
            clusterNum = nextid(model) # get a valid id for the model
            newCluster = Cluster(nextid(model), (0,0), [liposome.id, neighbor.id])
            add_agent_pos!(newCluster, model) # adds the cluster at the origin
                            
            # Make changes to the liposome and its neighbor
            liposome.type = :Member
            liposome.mems = [id]
            liposome.info = [newCluster.id, 2]
            neighbor.type = :Member
            neighbor.mems = [liposome.id]
            neighbor.info = [newCluster.id, 2]
            
            adjustCM(newCluster, model) # Moves the cluster to the correct position between the two liposomes
            break
        end
    end
    
    # Only do a random walk if the liposome has not bounded itself to a cluster.
    # If it has, then it's movement will be taken care of in cluster step. 
    if liposome.type == :Free
        liposome.vel = adjust_vel(liposome, model)
        displacement = liposome.vel .* model.dt
        walk!(liposome, displacement, model)
        #move_agent!(liposome, model, model.dt)
    end
end

function cluster_step!(cluster, model)
    # Cluster step first determines the cluster movement and rotation. Then it checks for any interaction between other
    # clusters, moves each member of the cluster, and then finally moves the cluster agent.
    
    cluster.vel = adjust_vel(cluster, model)
    # Movement of each member in the cluster
    for i in cluster.mems
        member = model[i]
                                        
        member.vel = cluster.vel
        displacement = member.vel .* model.dt
        walk!(member, displacement, model)
    end
                                    
    # Finally, move the cluster agent
    displacement = cluster.vel .* model.dt
    walk!(cluster, displacement, model)
                                    
    # If the cluster had no members, then it was absorbed earlier in system_step, so it should be removed
    if cluster.mems == []
        kill_agent!(cluster, model)
    end
end
                                
function cluster_interaction!(cluster, model)
    # Cluster step first determines the cluster movement and rotation. Then it checks for any interaction between other
    # clusters, moves each member of the cluster, and then finally moves the cluster agent.
    
    removed_clusters = []
                                
    # Loop through each member to see if any interact with any liposomes not part of the cluster
    for i in cluster.mems            
        member = model[i]
        neighbor_ids = nearby_ids(member, model, 1)
        
        for id in neighbor_ids
            # Check to see if this interaction results in aggregation
            if model[id].type != :Cluster
                if !successful_interaction(i, id, model)
                    continue # if not, skip this neighbor
                end
            else
                continue
            end
                                        
            neighbor = model[id]
                                        
            # Only have to check member-member interactions between members of different clusters because 
            # member-liposome interactions were already done first in the free liposome step.
            if neighbor.type == :Member && !(neighbor.id in member.mems)
                neighbor.mems = [neighbor.mems; member.id]
                member.mems = [member.mems; neighbor.id]
                                                
                if Int(neighbor.info[1]) != cluster.id
                    otherCluster = model[Int(neighbor.info[1])]
                    cluster.mems = [cluster.mems; otherCluster.mems] # Add the other cluster's members to this cluster's list

                    # Adjust the info of the new members to match this cluster                         
                    for i in otherCluster.mems
                        newMember = model[i]
                        newMember.info[1] = cluster.id
                    end
                    otherCluster.mems = [] # Set the other cluster's list to an empty vector
                    otherCluster.info = []  
                    push!(removed_clusters, otherCluster.id)
                    adjustCM(cluster, model) # Does final adjustments to the members and cluster agent's position
                end
            end
        end
    end
                                    
    calc_gyration_radius(cluster, model)
                                    
    if removed_clusters == []
        return -1
    else
        return removed_clusters
    end
end
                                

function cluster_rotate!(cluster, model)
    # Calculates a uniformly distributed random axis unit vector
    ϕ = 2*sqrt(cluster.info[3]*model.dt)*erfinv(2*rand(model.rng) - 1)
    if (rand() > 0.5)
        ϕ *= -1
    end
    rot = RotMatrix{2}(ϕ)
                                    
    # Movement of each member in the cluster
    for i in cluster.mems
        member = model[i]
                                        
        # First get the position of the member with respect to the cluster center of mass as bodyPos
        bodyPos = Vector{Float64}(collect(get_direction(cluster.pos, member.pos, model)))
        # Apply the quaternion rotation to bodyPos to find where that rotation would move that member to with respect to the cluster center of mass.
        newBodyPos = rot * bodyPos
        liposome_rotate!(i, model, ϕ)
        # The difference between the new bodyPos and the original bodyPos is how much the member is
        # shifted in any reference frame, not just with respect to the center of mass.
        shift = newBodyPos .- bodyPos
        walk!(member, Tuple(x for x in shift), model)
    end
end

function liposome_rotate!(id, model)
    ϕ = 2*sqrt(cluster.info[3]*model.dt)*erfinv(2*rand(model.rng) - 1)
    if (rand() > 0.5)
        ϕ *= -1
    end
    liposome_rotate!(id, model, ϕ)
end
                                
function liposome_rotate!(id, model, ϕ)
    model[id].oab[1] += ϕ
end
                                              
function system_step!(model)
    # This is the main step function for the program
    # Every step, first check for any interactions between liposomes
    # Then, move all liposomes and rotate all clusters
    # Finally, determine if any cluster moves 
                                    
    # Loop through all ids and make a vector of all free liposome and cluster agent ids
    freeIds = Vector{Int}()
    clusterIds = Vector{Int}()
    for id in allids(model)
        if model[id].type == :Free
            append!(freeIds, id)
        end
        if model[id].type == :Cluster
            append!(clusterIds, id)
        end
    end
                                            
    # Do the free liposome behavior first and then do cluster behavior
    for id in freeIds
        free_liposome_step!(model[id], model)
    end
                                            
    # Do cluster interactions
    for id in clusterIds
        removed_clusters = cluster_interaction!(model[id], model)
        # If clusters were absorbed in the step, remove it from the model now
        if removed_clusters != -1
            for remove in removed_clusters
                kill_agent!(model[remove], model)
                iter = 1
                for search in clusterIds
                    if search == remove
                        deleteat!(clusterIds, iter)
                        iter -= 1
                    end
                    iter+=1
                end
            end
        end
    end
    
    for id in clusterIds
        cluster_step!(model[id], model)
    end
    
    for id in clusterIds
        cluster_rotate!(model[id], model)
    end
                                            
    # Clear out any remaining dead clusters to avoid data collection errors
    for id in clusterIds
        if model[id].mems == []
            kill_agent!(model[id], model)
        end
    end
end

function adjustCM(cluster, model)
    # Function moves the cluster agent to the center of mass of the cluster.
    # Because of the boundaries being periodic, the best solution is to choose a cluster member's position
    # as the origin and then find the CM with respect to this origin, avoiding the boundary conditions.
    # This will likely fail if the cluster is longer than model.extent in any cardinal direction (which
    # should not happen)
    
    # Accumulator for the distance from the origin to all the members
    posCM = (0,0)
    # Use the first member as the origin of a different coordinate system
    origin = model[cluster.mems[1]].pos
    
    numMembers = size(cluster.mems)[1]
                                
    # Loops through all the cluster members
    for i in cluster.mems
        member = model[i]
        # Get the distance to the member from the origin, respecting periodicity
        posCM = posCM .+ get_direction(origin, member.pos, model)
        # Adjust each member's info to have the correct number of cluster members
        member.info[2] = numMembers
    end
    
    # Get CM in terms of our coordinate system
    posCM = posCM ./ numMembers
    
    # Move the cluster agent to this new CM
    move_agent!(cluster, origin, model)
    walk!(cluster, posCM, model)
    # Adjust statistics
    calc_gyration_radius(cluster, model)
    calc_hydrodynamic_radius(cluster, model)
    cluster.info[3] = k_b * model.T/ (8 * π * model.eta * (cluster.info[1] * 1e-7)^3 ) 
end

function successful_interaction(id_a, id_b, model)
    line = (model[id_b].pos .- model[id_a].pos)./2
    if abs(line[1]) > 1
        return false
    end
    angle_a = acos(line[1])
    angle_b = acos(-line[1])
    
    return in_mRNA(angle_a, id_a, model) ⊻ in_mRNA(angle_a, id_b, model)
end

function in_mRNA(ang, id, model)
    lip = model[id]
    if ang > ((lip.oab[1] + lip.oab[2]) % (2 * pi))
        return ang < ((lip.oab[1] + lip.oab[3]) % (2 * pi))
    else
        return ang > ((lip.oab[1] + lip.oab[3]) % (2 * pi))
    end
end
        
function calc_gyration_radius(cluster, model)
    gyration_ij = fill(length(cluster.mems)*2/5*(1e-7^2), (2,2))
    for i in cluster.mems
        r = get_direction(cluster.pos, model[i].pos, model) .* 1e-7
        r = vec([r[j] for j in 1:2])
        gyration_ij = gyration_ij + (dot(r,r) * I - reshape(kron(r,r), (2,2)))
    end
    cluster.info[1] = sqrt(sum(eigvals(gyration_ij)) ./ (1e-7^2))
end

function calc_hydrodynamic_radius(cluster, model)
    # Takes in a cluster and calculates the hydrodynamic radius for it, saving it in the appropriate location in cluster.info
    reciprocal_dist = 0.
    for i in 1:length(cluster.mems)
        for j in 1:(i-1)
            # Get the distance to each member from each other and add the reciprocal to a running sum, respecting periodicity
            reciprocal_dist += 1 / (euclidean_distance(model[cluster.mems[i]], model[cluster.mems[j]], model))
        end
    end
    r_hyd_reciprocal = reciprocal_dist / (2 * length(cluster.mems)^2)

    cluster.info[2] = 1 / r_hyd_reciprocal
end

function adjust_vel(agent, model)
    radius = 1e-7 # meters
    mass = model.mass # kg
    if agent.type == :Cluster
        radius = agent.info[1] * 1e-7
        mass *= length(agent.mems)
    end
        
    diff_const = k_b * model.T / (6 * π * model.eta * radius)
    lambda = mass * diff_const / (k_b * model.T * model.dt)
        
    new_vel = ((lambda/(1+lambda)) .* agent.vel) .+ ((1/(1+lambda)) .*
              (sqrt(2 * diff_const / model.dt) .* get_noise(model)) .+ 
              ((diff_const * mass / (k_b * model.T)) .* (0.,-9.81)))
    new_vel = new_vel ./ (1e-7) # convert back to units of 50 nm 
    return new_vel
end