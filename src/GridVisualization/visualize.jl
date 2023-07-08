using GraphRecipes
using Plots
using Graphs

function Plot_PowerGrid(grid ::PowerGrid; from_results=false, node_size=0.3,font_size=3)
    
    busbar_sections = []
    aux_buses = []
    reconf_lines = []
    couplers = []
    section2substation = Dict()
    for substation_id in keys(grid.Substations)
        if ! from_results
            non_aux_bus_in_substation = grid.Substations[substation_id].BusbarSections_IDs
            aux_buses_in_substation = grid.Substations[substation_id].Aux_Buses_IDs
            push!(section2substation,non_aux_bus_in_substation[1]=>substation_id)
            push!(section2substation,non_aux_bus_in_substation[2]=>substation_id)
            append!(busbar_sections,non_aux_bus_in_substation)
            append!(aux_buses,aux_buses_in_substation)
            append!(reconf_lines,grid.Substations[substation_id].Reconf_AuxLines_IDs)
            append!(couplers,grid.Substations[substation_id].Reconf_CouplerLines_IDs)
        else
            if grid.Substations[substation_id].is_split
                non_aux_bus_in_substation = grid.Substations[substation_id].BusbarSections_IDs
                aux_buses_in_substation = grid.Substations[substation_id].Aux_Buses_IDs
                push!(section2substation,non_aux_bus_in_substation[1]=>substation_id)
                push!(section2substation,non_aux_bus_in_substation[2]=>substation_id)
                append!(busbar_sections,non_aux_bus_in_substation)
                append!(aux_buses,aux_buses_in_substation)
                append!(reconf_lines,grid.Substations[substation_id].Reconf_AuxLines_IDs)
                append!(couplers,grid.Substations[substation_id].Reconf_CouplerLines_IDs)
            end
        end
    end

    names = Dict()
    nodeshapes = []
    nodecolors = []

    for bus in sort(collect(keys(grid.Buses)))
        if bus in busbar_sections
            if ! from_results
                if bus ∈ keys(grid.Substations)
                    push!(names,bus => string(section2substation[bus])*"-a")
                else
                    push!(names,bus => string(section2substation[bus])*"-b")
                end
                push!(nodeshapes, :hexagon)
                push!(nodecolors, colorant"#2986cc")
            else
                if grid.Substations[section2substation[bus]].is_split
                    if bus ∈ keys(grid.Substations)
                        push!(names,bus => string(section2substation[bus])*"-a")
                    else
                        push!(names,bus => string(section2substation[bus])*"-b")
                    end
                    push!(nodeshapes, :hexagon)
                    push!(nodecolors, colorant"#2986cc")
                else
                    push!(names,bus => string(bus)*" AC")
                    push!(nodeshapes, :circle)
                    push!(nodecolors, :lightblue)
                end
            end
        elseif bus in aux_buses
            push!(names,bus => string(bus)*" aux")
            push!(nodeshapes, :rect)
            push!(nodecolors, colorant"#2986cc")
        else
            push!(names,bus => string(bus)*" AC")
            push!(nodeshapes, :circle)
            push!(nodecolors, :lightblue)
        end
    end

    for bus in sort(collect(keys(grid.DCBuses)))
        push!(names,bus+grid.N_bus => string(bus)*" DC")
        push!(nodeshapes, :circle)
        push!(nodecolors, :lightblue)
    end

    if from_results
        g = Graphs.SimpleGraphs.SimpleGraph()
        add_vertices!(g,grid.N_bus+grid.N_dc_bus)
        edgestyles = Dict()
        edgecolors = Dict()

        for branch_id in keys(grid.Branches)
            my_branch = grid.Branches[branch_id]
            fr_bus_id = my_branch.Fr_bus_ID
            to_bus_id = my_branch.To_bus_ID
            add_edge!(g,fr_bus_id,to_bus_id)
            if branch_id ∈ reconf_lines
                
                if grid.Branches[branch_id].GeneralSwitch.SwitchingStatus == 0
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :dash)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :dash)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
                else
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :black)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :black)
                end
                
            elseif branch_id ∈ couplers
                if grid.Branches[branch_id].GeneralSwitch.SwitchingStatus == 0
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :dashdot)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :dashdot)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
                else
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :black)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :black)
                end
            else
                if grid.Branches[branch_id].GeneralSwitch.SwitchingStatus == 0
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
                else
                    push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                    push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                    push!(edgecolors, (fr_bus_id,to_bus_id) => :black)
                    push!(edgecolors, (to_bus_id,fr_bus_id) => :black)
                end
            end
        end

        for branch_id in keys(grid.DCBranches)
            my_branch = grid.DCBranches[branch_id]
            fr_bus_id = my_branch.Fr_bus_ID + grid.N_bus
            to_bus_id = my_branch.To_bus_ID + grid.N_bus
            add_edge!(g,fr_bus_id,to_bus_id)
            if grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus == 0
                push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
            else
                push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :blue)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :blue)
            end
        end

        for branch_id in keys(grid.Converters)
            my_branch = grid.Converters[branch_id]
            if my_branch.type == :ACDC
                fr_bus_id = my_branch.DC_Bus_ID + grid.N_bus
                to_bus_id = my_branch.AC_Bus_ID
            else
                fr_bus_id = my_branch.DC_Bus_ID
                to_bus_id = my_branch.AC_Bus_ID
            end

            add_edge!(g,fr_bus_id,to_bus_id)
            push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
            push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
            push!(edgecolors, (fr_bus_id,to_bus_id) => colorant"#9558B2")
            push!(edgecolors, (to_bus_id,fr_bus_id) => colorant"#9558B2")
        end

        for branch_id in keys(grid.DCLinks)
            my_branch = grid.DCLinks[branch_id]

            fr_bus_id = my_branch.Fr_bus_ID
            to_bus_id = my_branch.To_bus_ID

            add_edge!(g,fr_bus_id,to_bus_id)
            push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
            push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
            push!(edgecolors, (fr_bus_id,to_bus_id) => colorant"#9558B2")
            push!(edgecolors, (to_bus_id,fr_bus_id) => colorant"#9558B2")
        end

        for substation_id in keys(grid.Substations)
            if grid.Substations[substation_id].is_split
                fbus = grid.Substations[substation_id].BusbarSections_IDs[1]
                tbus = grid.Substations[substation_id].BusbarSections_IDs[2]
                add_edge!(g,fbus,tbus)
                push!(edgestyles,(fbus,tbus) => :dashdot)
                push!(edgestyles,(tbus,fbus) => :dashdot)
                push!(edgecolors, (fbus,tbus) => :grey)
                push!(edgecolors, (tbus,fbus) => :grey)
            end
        end


    else
        g = Graphs.SimpleGraphs.SimpleGraph()
        add_vertices!(g,grid.N_bus+grid.N_dc_bus)
        edgestyles = Dict()
        edgecolors = Dict()

        for branch_id in keys(grid.Branches)
            my_branch = grid.Branches[branch_id]
            fr_bus_id = my_branch.Fr_bus_ID
            to_bus_id = my_branch.To_bus_ID
            add_edge!(g,fr_bus_id,to_bus_id)
            if branch_id ∈ reconf_lines
                push!(edgestyles,(fr_bus_id,to_bus_id) => :dash)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :dash)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :black)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :black)
            elseif branch_id ∈ couplers
                push!(edgestyles,(fr_bus_id,to_bus_id) => :dashdot)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :dashdot)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
            else
                push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :black)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :black)
            end
        end

        for branch_id in keys(grid.DCBranches)
            my_branch = grid.DCBranches[branch_id]
            fr_bus_id = my_branch.Fr_bus_ID + grid.N_bus
            to_bus_id = my_branch.To_bus_ID + grid.N_bus
            add_edge!(g,fr_bus_id,to_bus_id)
            if grid.DCBranches[branch_id].GeneralSwitch.SwitchingStatus == 0
                push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :lightgray)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :lightgray)
            else
                push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
                push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
                push!(edgecolors, (fr_bus_id,to_bus_id) => :blue)
                push!(edgecolors, (to_bus_id,fr_bus_id) => :blue)
            end
        end

        for branch_id in keys(grid.Converters)
            my_branch = grid.Converters[branch_id]
            if my_branch.type == :ACDC
                fr_bus_id = my_branch.DC_Bus_ID + grid.N_bus
                to_bus_id = my_branch.AC_Bus_ID
            else
                fr_bus_id = my_branch.DC_Bus_ID
                to_bus_id = my_branch.AC_Bus_ID
            end

            add_edge!(g,fr_bus_id,to_bus_id)
            push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
            push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
            push!(edgecolors, (fr_bus_id,to_bus_id) => colorant"#9558B2")
            push!(edgecolors, (to_bus_id,fr_bus_id) => colorant"#9558B2")
        end

        for branch_id in keys(grid.DCLinks)
            my_branch = grid.DCLinks[branch_id]

            fr_bus_id = my_branch.Fr_bus_ID
            to_bus_id = my_branch.To_bus_ID

            add_edge!(g,fr_bus_id,to_bus_id)
            push!(edgestyles,(fr_bus_id,to_bus_id) => :solid)
            push!(edgestyles,(to_bus_id,fr_bus_id) => :solid)
            push!(edgecolors, (fr_bus_id,to_bus_id) => colorant"#9558B2")
            push!(edgecolors, (to_bus_id,fr_bus_id) => colorant"#9558B2")
        end

    end

    graphplot(g,names = names,
          nodeshape=nodeshapes, nodesize=node_size,
          fontsize=font_size,
          curves=false,
          color=:black,
          nodecolor=nodecolors,
          linewidth=2, edgestyle=edgestyles,
          edgecolor = edgecolors,method=:stress,trim=true)

    # return g,names,nodeshapes,nodecolors,edgestyles,edgecolors
end
