using JuMP, Gurobi

function solve_Optimal_TSPd(TT, DD, flying_range, bigM)
    model = Model(() -> Gurobi.Optimizer(gurobi_env))
    set_silent(model)
    m = size(TT)[1]
    n_nodes = m-2
    
    @variable(model, xT[1:m,1:m]>=0, Bin)
    @variable(model, xD[1:m,1:m]>=0, Bin)
    @variable(model, yT[2:m-1] >= 0, Bin) 
    @variable(model, yD[2:m-1] >= 0, Bin) 
    @variable(model, yC[2:m-1] >= 0, Bin) 
    @variable(model, a[1:m] >= 0) 

    @objective(model, Min, a[m])

    @constraint(model, sum(xT[1,j] for j=1:m ) == 1)   #(1d)
    @constraint(model, sum(xT[j,m] for j=1:m ) == 1)   #(1d)
    @constraint(model, sum(xD[1,j] for j=1:m ) == 1)   #(1g)
    @constraint(model, sum(xD[j,m] for j=1:m ) == 1)   #(1g)

    for i=2:m-1   #customer nodes
        @constraint(model, sum(xT[i,j] for j=1:m ) == sum(xT[j,i] for j=1:m ))   #(1b)
        @constraint(model, sum(xT[i,j] for j=1:m ) == yT[i]+yC[i] )              #(1c)
        @constraint(model, sum(xD[i,j] for j=1:m ) == sum(xD[j,i] for j=1:m ))   #(1e)
        @constraint(model, sum(xD[i,j] for j=1:m ) == yD[i]+yC[i] )              #(1f)
        @constraint(model, yT[i] + yD[i] + yC[i] == 1)                           #(1h)
        @constraint(model, xD[1,i] + xD[i,m] <= 1)                               #(1l)
        for j=i+1:m-1
            @constraint(model, xD[i,j] + xD[j,i] <= yC[i] + yC[j])               #(1k)
        end 
        @constraint(model, sum(DD[j,i]*xD[j,i] for j=1:m) + sum(DD[i,k]*xD[i,k] for k=1:m) <= flying_range + bigM*(1-yD[i]))   #(5d)
    end
    # for j=2:m-1
    #     @constraint(model, sum(xT[i,j] for i=2:m-1 ) == 1)
    # end

    for i=1:m   # all nodes (depot + customers + depot)
        for j=1:m
            @constraint(model, a[i] + TT[i,j] <= a[j] + bigM*(1-xT[i,j]))    #(1i)
            @constraint(model,a[i] + DD[i,j] <= a[j] + bigM*(1-xD[i,j]))     #(1j)
            if DD[i,j] > flying_range
                @constraint(model, xD[i,j] <= xT[i,j])                              #(5a)
            end
            if i == j
                @constraint(model, xT[i,j] == 0)
                @constraint(model, xD[i,j] == 0)
            end

        end
    end

    @constraint(model, sum(TT[i,j]*xT[i,j] for i=1:m, j=1:m)<= a[m])     #(Valid enqualities)
    @constraint(model, sum(DD[i,j]*xD[i,j] for i=1:m, j=1:m)<= a[m])     #(Valid enqualities)
#     set_time_limit_sec(model, 300)
    optimize!(model)


    #Initialize variables and save outputs
    xT_out = zeros(m,m)
    xD_out = zeros(m,m)
    yT_out = zeros(n_nodes)    
    yD_out = zeros(n_nodes) 
    yC_out = zeros(n_nodes) 
    a_out = zeros(m)
    for i=1:m
        a_out = value(a[i])
        for j=1:m
            xT_out[i,j]= value(xT[i,j])
            xD_out[i,j]= value(xD[i,j])
        end
    end

    for i=1:m-2

        yT_out[i]= value(yT[i+1])
        yD_out[i]= value(yD[i+1])
        yC_out[i]= value(yC[i+1])

    end

    obj_value = objective_value(model)
    return xT_out, xD_out, yT_out, yD_out, yC_out, a_out, obj_value
end
