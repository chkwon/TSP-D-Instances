using TSPDrone 

# filenames = ["uniform-$(60+i)-n20.txt" for i in 1:1]
# filenames = ["uniform-$i-n11.txt" for i in 1:10]
filenames = ["uniform-$(50+i)-n10.txt" for i in 1:10]

function test(filenames)
    objs = Float64[]
    for file in filenames
        x, y, truck_cost_factor, drone_cost_factor = read_data("uniform", file)

        objective_value, truck_route, drone_route = solve_tspd(x, y, truck_cost_factor, drone_cost_factor; )

        @show objective_value, truck_route, drone_route
        push!(objs, objective_value)
    end

    mean(objs)
end

test(filenames)