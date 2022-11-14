using TSPDrone 
using TSPLIB
using Concorde
using Statistics

function print_time(total_sec)
    h = floor(Int, total_sec / 3600)
    m = floor(Int, (total_sec - h * 3600) / 60)
    s = floor(Int, (total_sec - h * 3600 - m * 60))
    ss = round(total_sec * 100) / 100

    if h > 0
        return "$(h)h $(m)m"
    elseif m > 0
        return "$(m)m $(s)s"
    elseif s > 0
        return "$(s)s"
    else
        return "$(ss)s"
    end
end

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

function test_agatz(dist, n_nodes; DPS25=true, alpha=2, n_samples=0, device="cpu")
    prefix = alpha == 2 ? "" : "-alpha_$(alpha)"
    
    filenames = Dict()
    filenames[10] = ["$(dist)$(prefix)-$(50 + i)-n10.txt" for i in 1:10]
    filenames[20] = ["$(dist)$(prefix)-$(60 + i)-n20.txt" for i in 1:10]
    filenames[50] = ["$(dist)$(prefix)-$(70 + i)-n50.txt" for i in 1:10]
    filenames[75] = ["$(dist)$(prefix)-$(80 + i)-n75.txt" for i in 1:10]
    filenames[100] = ["$(dist)$(prefix)-$(90 + i)-n100.txt" for i in 1:10]
    filenames[175] = ["$(dist)$(prefix)-$(100 + i)-n175.txt" for i in 1:10]
    filenames[250] = ["$(dist)$(prefix)-$(110 + i)-n250.txt" for i in 1:10]

    for n in n_nodes
        if DPS25
            n_grp = n < 25 ? 1 : Int(n / 25)
        else 
            n_grp = 1
        end
        
        objs = Float64[]
        objs_RL = Float64[]

        if n_samples == 0
            t = time()
            for filename in filenames[n]
                x, y, truck_cost_factor, drone_cost_factor = read_data(dist, filename)
                @assert truck_cost_factor / drone_cost_factor == alpha

                # manhattan_dist_mtx = [abs(x[i]-x[j]) + abs(y[i]-y[j]) for i in 1:n, j in 1:n]
                # euclidean_dist_mtx = [sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) for i in 1:n, j in 1:n]
                
                # truck_dist_mtx = euclidean_dist_mtx ./ 40 .* 100
                # drone_dist_mtx = manhattan_dist_mtx ./ 40 ./ 2 .* 100

                # objective_value, truck_route, drone_route = solve_tspd(truck_dist_mtx, drone_dist_mtx; n_groups=n_grp, method="TSP-ep-all")

                result = solve_tspd(x, y, truck_cost_factor, drone_cost_factor; n_groups=n_grp, method="TSP-ep-all")
                obj = result.total_cost
                et = print_time(time() - t)
                @show filename, obj, n_grp, et
                push!(objs, obj)
            end
            t_dps = time() - t

            open("$(dist)$(prefix)-n$n-DPS-n_groups_$(n_grp).txt", "w") do io 
                println(io, objs)
                println(io, mean(objs))
                println(io, t_dps / 10)
            end
        end

        if n_samples > 0

            t = time()
            for filename in filenames[n]
                x, y, truck_cost_factor, drone_cost_factor = read_data(dist, filename)

                @assert truck_cost_factor / drone_cost_factor == 2.0
                obj_value, _, _ = solve_tspd_RL(x, y, n_samples=n_samples, device=device)
                @show filename, obj_value, n_grp
                push!(objs_RL, obj_value[1])
            end
            t_rl = time() - t

            open("$(dist)$(prefix)-n$n-RL-n_samples_$(n_samples).txt", "w") do io
                println(io, objs_RL)
                println(io, mean(objs_RL))
                println(io, t_rl / 10)
            end

        end

        println("n = $(n), mean(DPS25)= $(mean(objs)), mean(RL) = $(mean(objs_RL))")

    end

end

# test_agatz(n_samples=1, device="cpu")



# alpha = 2
#test_agatz("uniform", [10, 20, 50, 75, 100, 175, 250]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("singlecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("doublecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")

# alpha = 2
# #test_agatz("uniform", [10, 20, 50, 75, 100, 175, 250]; alpha=alpha, n_samples=0, device="cuda")
# # test_agatz("singlecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# # test_agatz("doublecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")

# alpha = 1
# # test_agatz("uniform", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("singlecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("doublecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")

alpha = 3
# test_agatz("uniform", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("singlecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("doublecenter", [10, 20, 50, 75, 100]; alpha=alpha, n_samples=0, device="cuda")
# test_agatz("doublecenter", [10, 100]; alpha=alpha, n_samples=0, device="cuda")

# test_agatz("uniform", [20, 50, 100]; n_samples=1, device="cuda")
# test_agatz("singlecenter", [20, 50, 100]; n_samples=1, device="cuda")
# test_agatz("doublecenter", [20, 50, 100]; n_samples=1, device="cuda")


# Please run these:
# test_agatz("uniform", [20, 50, 100]; n_samples=1, device="cuda")
# test_agatz("uniform", [20, 50, 100]; n_samples=4800, device="cuda")


# Table 4
alpha = 2
# test_agatz("uniform", [20, 50, 100]; DPS25=false, alpha=alpha, n_samples=0, device="cpu")
# test_agatz("uniform", [20, 50, 100]; DPS25=true, alpha=alpha, n_samples=0, device="cpu")
# test_agatz("uniform", [20, 50, 100]; n_samples=1, device="cpu")
# test_agatz("uniform", [20, 50, 100]; n_samples=4800, device="cpu")
test_agatz("uniform", [100]; n_samples=4800, device="cuda")