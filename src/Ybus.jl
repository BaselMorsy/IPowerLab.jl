function Y_Bus(LineData,N)

    Y_bus = zeros(N,N)*0*im
    M = size(LineData,1);
    otherBus = 0;
    for i in 1:N
        BusIndex = i;
        for j in 1:M
            b1 = LineData[j,1];
            b2= LineData[j,2];
            b = [b1 b2];

            if b1 == BusIndex || b2 == BusIndex
                for k in 1:2
                    if b[k] != BusIndex
                        otherBus = b[k];
                    end
                end

                r = LineData[j,3];
                x = LineData[j,4];

                Z = r + x*im;

                Y_sh = LineData[j,5]*im;
                Y_bus[BusIndex,BusIndex] = Y_bus[BusIndex,BusIndex] + 1/Z + 0.5*Y_sh;
                Y_bus[BusIndex,otherBus] = Y_bus[BusIndex,otherBus] - 1/Z
            end
        end
    end

    b = zeros(N,N)
    for l in 1:M
        b[LineData[l,:fbus],LineData[l,:tbus]] = LineData[l,:b]
        b[LineData[l,:tbus],LineData[l,:fbus]] = LineData[l,:b]
    end

     return Y_bus, b
 end


 function Y_Bus_Grid(grid ::PowerGrid)

    @assert length(keys(grid.Buses)) == grid.N_bus

    branch_dictionary = Dict()
    for branch_id in keys(grid.Branches)
        push!(branch_dictionary, Set([grid.Branches[branch_id].Fr_bus_ID,grid.Branches[branch_id].To_bus_ID]) => branch_id)
    end

    y_bus = zeros(grid.N_bus,grid.N_bus)*0*im
    b_line = zeros(grid.N_bus,grid.N_bus)*0*im
    for i in keys(grid.Buses)
        for j in keys(grid.Buses)
            if i != j
                if Set([i,j]) in keys(branch_dictionary)
                    my_branch = grid.Branches[branch_dictionary[Set([i,j])]]
                    r = my_branch.r
                    x = my_branch.x
                    z = r + x*im
                    b = my_branch.b
                    b_line[i,i] += im*b/2
                    y_bus[i,j] = -1/z
                end
            end
        end
    end
    
    for i in keys(grid.Buses)
        
        sums = collect(setdiff(Set(keys(grid.Buses)),Set([i])))
        y_bus[i,i] = b_line[i,i] + sum(y_bus[i,k] for k in sums)
    end

    return y_bus,b_line
 end