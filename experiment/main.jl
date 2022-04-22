using TSPDrone 
using TSPLIB
using Concorde

function read_data(subfolder, instance)
    file = open(joinpath(@__DIR__, "..", subfolder, instance), "r")

    readline(file) # /*The speed of the Truck*/
    truck_cost_factor = parse(Float64, readline(file))
    
    readline(file) # /*The speed of the Drone*/
    drone_cost_factor = parse(Float64, readline(file))
    
    readline(file) # /*Number of Nodes*/
    n_nodes = parse(Int, readline(file))
    
    readline(file) # /*The Depot*/
    depot = split(readline(file), " ")
    depot_coordinates = parse.(Float64, depot[1:2])
    
    readline(file) # /*The Locations (x_coor y_coor name)*/
    customer_coordinates = Matrix{Float64}(undef, n_nodes-1, 2)
    for i in 1:n_nodes-1
        customer = split(readline(file), " ")
        customer_coordinates[i, :] = parse.(Float64, customer[1:2])   
    end
    
    x = vcat(depot_coordinates[1], customer_coordinates[:, 1])
    y = vcat(depot_coordinates[2], customer_coordinates[:, 2])

    close(file)

    return x, y, truck_cost_factor, drone_cost_factor
end

function test_TSPLIB100()
    instances = [:kroA100, :kroC100, :kroD100, :kroE100, :rd100]
    for inst in instances
        tsp = readTSPLIB(inst)
        x = tsp.nodes[:, 1]
        y = tsp.nodes[:, 2]
        truck_cost_factor = 1.0
        drone_cost_factor = 0.5

        manhattan_dist_mtx = [abs(x[i]-x[j]) + abs(y[i]-y[j]) for i in 1:tsp.dimension, j in 1:tsp.dimension]
        dist_mtx = round.(Int, manhattan_dist_mtx ./ 40 .* 100)
        tsp_tour, tsp_len = solve_tsp(dist_mtx)

        drone_dist_mtx = tsp.weights ./40 .* 100
        truck_dist_mtx = manhattan_dist_mtx ./ 40 .* 100

        flying_range = 40 * 40/60 * 100

        objective_value, truck_route, drone_route = solve_tspd(truck_dist_mtx, drone_dist_mtx, flying_range=flying_range, n_groups=10, method="TSP-ep-all")

        println("$(inst), tsp=$(tsp_len), tspd=$(objective_value)")

    end
end

function test_agatz(;n_samples=1, device="cpu")
    filenames = Dict()
    filenames[20] = ["uniform-$(60 + i)-n20.txt" for i in 1:10]
    filenames[50] = ["uniform-$(70 + i)-n50.txt" for i in 1:10]
    filenames[100] = ["uniform-$(90 + i)-n100.txt" for i in 1:10]

    n_groups = Dict()
    n_groups[20] = 1
    n_groups[50] = 2
    n_groups[100] = 4

    n_nodes = [20, 50, 100]

    for n in n_nodes
        n_grp = n_groups[n]
        objs = Float64[]
        objs_RL = Float64[]

        t = time()

        for filename in filenames[n]
            @info("Solving $filename by DPS")
            x, y, truck_cost_factor, drone_cost_factor = read_data("uniform", filename)
            
            objective_value, truck_route, drone_route = solve_tspd(x, y, truck_cost_factor, drone_cost_factor, n_groups=n_grp, method="TSP-ep-all")
            push!(objs, objective_value)
        end
        t_dps = time() - t

        open("uniform-n$n-DPS.txt", "w") do io 
            println(io, objs)
            println(io, mean(objs))
            println(io, t_dps / 10)
        end


        t = time()
        for filename in filenames[n]
            @info("Solving $filename by RL")
            x, y, truck_cost_factor, drone_cost_factor = read_data("uniform", filename)

            @assert truck_cost_factor / drone_cost_factor == 2.0
            obj_value, truck_route, drone_route = solve_tspd_RL(x, y, n_samples=n_samples, device=device)
            push!(objs_RL, obj_value[1])
        end
        t_rl = time() - t

        open("uniform-n$n-RL.txt", "w") do io
            println(io, objs_RL)
            println(io, mean(objs_RL))
            println(io, t_rl / 10)
        end

        println("n = $(n), mean(DPS)= $(mean(objs)), mean(RL) = $(mean(objs_RL))")

    end

end

test_agatz(n_samples=1, device="cpu")


# Please run 
# test_agatz(n_samples=1, device="cpu")
