#Create arc set
function create_arc(V)
    N = length(V)
    #complete graph
    A = Set()
    for i in 1:N
        for j in 1:N
            if i!=j
                push!(A, (i,j))
            end
        end
    end
    #Add depot node i in V\{1} -> N+1
    for i in 2:N
        push!(A, (i,N+1))
    end
    return A
end

#return TSP solution from route variable
function sol_TSP(A, X, N) 
    p = Set()
    for (i,j) in A
        if abs(X[(i,j)]-1.0) < 1e-5
            if j == N+1
                push!(p, (i, N+1) )
            else 
                push!(p, (i,j))
            end
        end
    end

    start = []
    end_v = []
    for (i,j) in p
        push!(start,i)
        push!(end_v,j)
    end

    path = [1]
    next = findfirst(isequal(1), start)
    push!(path, end_v[next])

    for i in 1:length(start)-1
        now = path[end]
        next = findfirst(isequal(now), start)
        push!(path, end_v[next])
    end
    return path
end

#solve it using formulation by Gurobi
using JuMP, Gurobi
function formulation(truck_cost_mtx, drone_cost_mtx, N, e, M)
    tT = hcat(truck_cost_mtx, truck_cost_mtx[:,1]) #add colummn i->N+1
    tD = hcat(drone_cost_mtx, drone_cost_mtx[:,1]) 
    #M = 1e3 #big-M
    
    V = Set(collect(1:N)) #V = {1,2,...,N}
    V1 = delete!(deepcopy(V),1) # V1 = {2,3,..,N}
    V2 = deepcopy(V)
    push!(V2, N+1) # V2={1,2,...,N+1}
    A = create_arc(V)

    
    m = Model(() -> Gurobi.Optimizer(gurobi_env))
    # set_optimizer_attribute(m, "TimeLimit", Time_limit)
    #set_optimizer_attribute(m, "IntFeasTol", 1e-9)
    set_optimizer_attribute(m, "Outputflag", 0)
    
    @variable(m, xT[(i,j) in A], Bin)
    @variable(m, xD[(i,j) in A], Bin)
    @variable(m, nT[i in V1]>=0, Int)
    @variable(m, nD[i in V1]>=0, Int)
    @variable(m, yT[i in V1], Bin)
    @variable(m, yD[i in V1], Bin)
    @variable(m, yC[i in V1], Bin)
    @variable(m, a[i in V2]>=0)
    @variable(m, α[i in V1,j in V1], Bin)
    @variable(m, f[i in V2]>=0)
    
    @objective(m, Min, a[N+1])
    
    #Flow in the Truck Route
    @constraint(m, sum(xT[(i,j)] for (i,j) in A if i==1) == 1) 
    @constraint(m, sum(xT[(i,j)] for (i,j) in A if j==N+1) == 1)
    @constraint(m, [i in V1], sum(xT[(ii,j)] for (ii,j) in A if ii==i) - sum(xT[(j,ii)] for (j,ii) in A if ii==i) == 0)
    
    #Flow in the Drone Route
    @constraint(m, sum(xD[(i,j)] for (i,j) in A if i==1) == 1) 
    @constraint(m, sum(xD[(i,j)] for (i,j) in A if j==N+1) == 1)
    @constraint(m, [i in V1], sum(xD[(ii,j)] for (ii,j) in A if ii==i) - sum(xD[(j,ii)] for (j,ii) in A if ii==i) == 0)
    
    #Node category
    @constraint(m, [i in V1], yT[i] + yD[i] + yC[i] == 1)
    for i in V1
        for j in V1
            if i!=j
                @constraint(m, xD[(i,j)] + xD[(j,i)] <= yC[i] + yC[j] + 2*(1-α[i,j]) )
                @constraint(m, xD[(i,j)] + xD[(j,i)] <= yD[i] + yD[j] + 2*α[i,j] ) 
            end
        end
    end
    @constraint(m, [i in V1], sum(xT[(ii,j)] for (ii,j) in A if ii==i) == yT[i] + yC[i])
    @constraint(m, [i in V1], sum(xD[(ii,j)] for (ii,j) in A if ii==i) == yD[i] + yC[i])
    
    #Subtour elimination
    for i in V1
        for j in V1
            if i!=j
            @constraint(m, nT[j] >= nT[i] + N*xT[(i,j)] - (N-1)) #subtour elimination for the truck route
            @constraint(m, nD[j] >= nD[i] + N*xD[(i,j)] - (N-1)) #subtour elimination for the drone route
            end
        end
    end
    
    #Arrival Time
    @constraint(m, [(i,j) in A], a[i] + tT[i,j] <= a[j] + M*(1-xT[(i,j)]))
    @constraint(m, [(i,j) in A], a[i] + tD[i,j] <= a[j] + M*(1-xD[(i,j)]) + M*xT[(i,j)])
    @constraint(m, sum(tT[i,j]*xT[(i,j)] for (i,j) in A) <= a[N+1])
    @constraint(m, sum(tD[i,j]*xD[(i,j)] for (i,j) in A) <= a[N+1])
    
    #Drone flying Range
    for (i,j) in A
        if tD[i,j] > e
             @constraint(m, xD[(i,j)] <= xT[(i,j)])
         end
         @constraint(m, f[j] >= f[i] + tD[i,j]*xD[(i,j)] - M*(1-xD[(i,j)]) - M*xT[(i,j)])
     end
 
     for i in V2
         @constraint(m, f[i] <= e)
     end



    time = @timed optimize!(m)
    #@show raw_status(m)
    Tr = value.(xT)
    Dr = value.(xD)
    yT = value.(yT)
    yD = value.(yD)
    yC = value.(yC)
    arr = value.(a)
    sol = objective_value(m)

    tr, dr = sol_TSP(A, Tr, N), sol_TSP(A, Dr, N) #route 

    return time, sol, tr, dr
end


# Gurobi_time, Gurobi_sol, Gurobi_tr, Gurobi_dr = formulation(truck_cost_mtx, drone_cost_mtx, N, e, Time_limit, M)
#V = {1,2,...,N}

# ####################Input########################
# N = 11         #number of nodes include the depot
# vt = 1.0       #truck speed
# vd = 2.0       #drone speed
# e = 100/vd     #drone flight time = range / vd
# ###################################################