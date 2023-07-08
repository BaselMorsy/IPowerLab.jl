using Parameters
using DataFrames
using CSV
using Graphs
using GraphRecipes
using Plots

include("Components/Grid.jl")
include("Components/Bus.jl")
include("Components/Branch.jl")
include("Components/Generator.jl")
include("Components/Load.jl")
include("Components/Substation.jl")
include("Components/Storage.jl")
include("Components/Switch.jl")
include("Components/DCBranch.jl")
include("Components/DCBus.jl")
include("Components/Converter.jl")
include("Components/DCLink.jl")
include("Ybus.jl")


function BuildGrid(Sbase,LineData,BusData,GenData)
    """
    Soon will be removed!
    """
    Y_bus = []
    b_line = []
    Buses = Dict()
    Branches = Dict()
    Loads = Dict()
    Generators = Dict()

    gens_at_bus = Dict()
    branches_at_bus = Dict()

    switch_counter = 0
    branch_counter = 0
    bus_counter = 0
    load_counter = 0
    gen_counter = 0

    N_bus = []
    N_gen = []
    N_load = []
    N_branch = []
    N_switch = []

    if LineData != []
        N_bus = length(unique!(vcat(unique!(Vector(LineData[!,:fbus])),unique!(Vector(LineData[!,:tbus])))))
        Y_bus,b_line = Y_Bus(LineData,N_bus)

        for line in eachrow(LineData)
            objtyp = 2
            branch_counter += 1
            switch_counter += 1
            new_switch = Switch(SwitchID=switch_counter,Switch_ObjID=(objtyp,line[:ID]))
            new_line = Branch(LineID=line[:ID],Fr_bus_ID=line[:fbus],To_bus_ID=line[:tbus],
                BranchType=0,GeneralSwitch=new_switch,r=line[:r],x=line[:x],b=line[:b],
                rating=line[:rate])
            push!(Branches,line[:ID]=>new_line)

            if new_line.Fr_bus_ID in keys(branches_at_bus)
                push!(branches_at_bus[new_line.Fr_bus_ID],new_line.LineID)
            else
                push!(branches_at_bus,new_line.Fr_bus_ID=>[new_line.LineID])
            end

            if new_line.To_bus_ID in keys(branches_at_bus)
                push!(branches_at_bus[new_line.To_bus_ID],new_line.LineID)
            else
                push!(branches_at_bus,new_line.To_bus_ID=>[new_line.LineID])
            end
        end

    end
    
    if GenData != []
        for gen in eachrow(GenData)
            objtyp = 3
            gen_counter += 1
            switch_counter += 1
            new_switch = Switch(SwitchID=switch_counter,Switch_ObjID=(objtyp,gen_counter))
            new_gen = Generator(GenID=gen_counter,GenBus_ID=gen[:bus],GeneralSwitch=new_switch,C0=gen[:C0],C1=gen[:C1],C2=gen[:C2],
                Pg_max=gen[:Pmax],Pg_min=gen[:Pmin],Qg_max=gen[:Qmax],Qg_min=gen[:Qmin])
            push!(Generators,gen_counter=>new_gen)

            if new_gen.GenBus_ID in keys(gens_at_bus)
                push!(gens_at_bus[new_gen.GenBus_ID],new_gen.GenID)
            else
                push!(gens_at_bus,new_gen.GenBus_ID=>[new_gen.GenID])
            end

        end

    end

    if BusData != []
        
        for bus in eachrow(BusData)
            switch_counter +=1
            bus_counter +=1
            load_counter +=1
            objtyp = 1
            new_switch = Switch(SwitchID=switch_counter,Switch_ObjID=(objtyp,bus[:bus_i]))
            new_load = Load(LoadID=bus[:bus_i],LoadBus_ID=bus[:bus_i],Pd=bus[:Pd],Qd=bus[:Qd],Shedding_Cost=1e4, Pd_t=[bus[:Pd]])
            new_bus = Bus(BusID=bus[:bus_i],GeneralSwitch=new_switch,V_max=bus[:Vmax],V_min=bus[:Vmin])
            
            if new_bus.BusID in keys(branches_at_bus)
                new_bus.ConnectedLinesIDs = branches_at_bus[new_bus.BusID]
            end

            new_bus.ConnectedLoadsIDs = [new_load.LoadID]

            if new_bus.BusID in keys(gens_at_bus)
                new_bus.ConnectedGensIDs = gens_at_bus[new_bus.BusID]
            end

            push!(Buses, bus[:bus_i]=>new_bus)
            push!(Loads,bus[:bus_i]=>new_load)
        end
    end

    N_gen = gen_counter
    N_load = load_counter
    N_branch = branch_counter
    N_switch = switch_counter

    grid = PowerGrid(GridID=1,Areas=1)
    grid.Buses = Buses
    grid.Loads = Loads
    grid.Branches = Branches
    grid.Generators = Generators

    grid.N_bus = N_bus
    grid.N_gen = N_gen
    grid.N_load = N_load
    grid.N_branch = N_branch
    grid.N_switch = N_switch

    grid.Y_bus =  Y_bus
    grid.b_line = b_line
    grid.S_base = Sbase

   return grid
end

function show_cases(verbose; type=:AC_cases)

    """
    type: could be :AC_cases for pglib_opf benchmark cases,
        or :ACDC_cases for hybrid benchmark systems, or :UC_cases for unit commitment cases
    """
    
    lst_dir = abspath(joinpath(@__DIR__, ".."))
    if type == :AC_cases
        lst = readdir(string(lst_dir,"\\test\\Systems\\pglib-opf\\benchmark_cases"))
    elseif type == :ACDC_cases
        lst = readdir(string(lst_dir,"\\test\\Systems\\ACDC"))
    elseif type == :UC_cases
        lst = readdir(string(lst_dir,"\\test\\Systems\\AC_UC"))
        if length(lst) == 0
            println("No UC cases downloaded, please use `download_UC_case(case_name::string; date::string)`")
        end
    end

    if verbose
        for i in 1:length(lst)
            println("[INFO] ", lst[i]," - ", i)
        end
    end
    return lst
end

function load_system(case_id; type=:AC_cases, date="2017-02-18")
    """
    date: only relevant for UC_cases. Use a date of a downloaded case
    """
    lst_dir = abspath(joinpath(@__DIR__, ".."))
    lst = show_cases(false, type=type)
    CASE_NAME= lst[case_id]

    if type == :AC_cases
        CASE_DIR = string(lst_dir,"\\test\\Systems\\pglib-opf\\benchmark_cases")
        m_file_path = joinpath(CASE_DIR, CASE_NAME)
        return parse_matpower_case(m_file_path; start_node_from_1=true)
    elseif type == :ACDC_cases
        CASE_DIR = string(lst_dir,"\\test\\Systems\\ACDC")
        m_file_path = joinpath(CASE_DIR, CASE_NAME)
        return parse_matpower_case(m_file_path; start_node_from_1=true)
    elseif type == :UC_cases
        CASE_DIR = string(lst_dir,"\\test\\Systems\\AC_UC")
        json_grid_file_path = joinpath(joinpath(CASE_DIR, CASE_NAME), date*".json.gz")
        std_case_lst = show_cases(false, type=:AC_cases)
        standard_case_name = []
        special_suffixes = ["pegase", "rte", "wp", "sp", "wop", "sop"]

        for suffix in special_suffixes
            if occursin(suffix, CASE_NAME)
                CASE_NAME = replace(CASE_NAME, suffix => "_$suffix")
                break
            end
        end
        i = 0
        found = true
        for ac_case in std_case_lst
            i += 1
            if occursin(CASE_NAME, ac_case)
                standard_case_name = (ac_case, i)
                found = true
                break
            else
                found = false
            end
        end

        if found
            return parse_json_case(json_grid_file_path,standard_case_name)
        else
            println("Could not find an equivalent system in pglib cases! Aborting!")
            return -1
        end
    end

end

function load_case(case_id,Sbase)
    """
    Soon will be removed!
    """
   lst = show_cases(false)
   CASE_DIR = string(string(@__DIR__,"\\Systems\\csv"),"/",lst[case_id])

   line_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","branch.csv"),DataFrame ;copycols = true))
   LineData = DataFrame(fbus = line_data_raw[:,1],tbus = line_data_raw[:,2],r = line_data_raw[:,3],x = line_data_raw[:,4],b = line_data_raw[:,5],rate = line_data_raw[:,6]./Sbase,ID = collect(1:length(line_data_raw[:,6])), is_aux = zeros(length(line_data_raw[:,1]),))

   bus_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","bus.csv"),DataFrame ;copycols = true))
   BusData = DataFrame(bus_i = bus_data_raw[:,1],type = bus_data_raw[:,2],Pd = bus_data_raw[:,3],Qd = bus_data_raw[:,4],Vmax = bus_data_raw[:,12],Vmin = bus_data_raw[:,13],is_aux = zeros(length(bus_data_raw[:,1]),))

   gen_raw =  DataFrame(CSV.read(string(CASE_DIR,"/","gen.csv"),DataFrame ;copycols = true))
   gen_cost_raw =  DataFrame(CSV.read(string(CASE_DIR,"/","gencost.csv"),DataFrame ;copycols = true))
   GenData = DataFrame(bus = gen_raw[:,1],C2 = gen_cost_raw[:,5],C1 = gen_cost_raw[:,6],C0 = gen_cost_raw[:,7],Pmax = gen_raw[:,9],Pmin = gen_raw[:,10], Qmax = gen_raw[:,4], Qmin = gen_raw[:,5])


   return BuildGrid(Sbase,LineData,BusData,GenData)
end

function load_saved_case(case_id)
    lst = show_cases(false)

    if occursin("pglib_opf_",lst[case_id])
        file_name = replace(lst[case_id],"pglib_opf_"=>"")
    elseif occursin("UC_",lst[case_id])
        file_name = replace(lst[case_id],"UC_"=>"")
    end
    CASE_DIR = string(string(@__DIR__,"\\Systems\\csv"),"/",lst[case_id],"/",file_name)
    return load_grid_object(CASE_DIR)
end

function set_load_multiplier!(grid ::PowerGrid,value)

    for key in keys(grid.Loads)
        grid.Loads[key].Pd *= value/grid.load_multiplier
        grid.Loads[key].Qd *= value/grid.load_multiplier
        if grid.Loads[key].Pd_t != []
            grid.Loads[key].Pd_t .*= value/grid.load_multiplier
        end
        if grid.Loads[key].Qd_t != []
            grid.Loads[key].Qd_t .*= value/grid.load_multiplier
        end
    end

    grid.load_multiplier = value
end

function set_load_profile_multiplier!(grid ::PowerGrid,value)
    
    for key in keys(grid.Loads)
        grid.Loads[key].Pd_t *= value/grid.load_multiplier
        grid.Loads[key].Qd_t *= value/grid.load_multiplier
    end

    grid.load_multiplier = value
end

function set_generation_multiplier!(grid ::PowerGrid,value)

    for key in keys(grid.Generators)
        grid.Generators[key].Pg_max *= value/grid.gen_multiplier
        grid.Generators[key].Qg_max *= value/grid.gen_multiplier
    end

    grid.gen_multiplier = value
end

function set_line_capacity_multiplier!(grid::PowerGrid, value)

    for key in keys(grid.Branches)
        grid.Branches[key].rating *= value / grid.line_capacity_multiplier
    end
    grid.line_capacity_multiplier = value
end

function set_converter_capacity_multiplier!(grid::PowerGrid, value)
    for key in keys(grid.Converters)
        grid.Converters[key].rate *= value / grid.converter_capacity_multiplier
        gen_dc_id = grid.Converters[key].gen_dc_id
        gen_ac_id = grid.Converters[key].gen_ac_id
        grid.Generators[gen_dc_id].Pg_max *= value / grid.converter_capacity_multiplier
        grid.Generators[gen_dc_id].Qg_max *= value / grid.converter_capacity_multiplier

        grid.Generators[gen_ac_id].Pg_max *= value / grid.converter_capacity_multiplier
        grid.Generators[gen_ac_id].Qg_max *= value / grid.converter_capacity_multiplier

        grid.Generators[gen_dc_id].Pg_min *= value / grid.converter_capacity_multiplier
        grid.Generators[gen_dc_id].Qg_min *= value / grid.converter_capacity_multiplier

        grid.Generators[gen_ac_id].Pg_min *= value / grid.converter_capacity_multiplier
        grid.Generators[gen_ac_id].Qg_min *= value / grid.converter_capacity_multiplier
    end

    

    grid.converter_capacity_multiplier = value
end

function convert_bus2substation!(grid ::PowerGrid,BusIDs=[],add_b2b_VSC=false,VSC_capacity=10)
    N_Sections=2
    for bus in BusIDs
        #= Don't forget to clear the main bus from loads, gens, lines and storage !!! =#
        Bus_obj = grid.Buses[bus]
        N_elements_connected = length(Bus_obj.ConnectedGensIDs) + length(Bus_obj.ConnectedLinesIDs) + length(Bus_obj.ConnectedLoadsIDs) + length(Bus_obj.ConnectedStorageIDs)
        N_additional_nonaux_buses = N_Sections - 1

        gen_ID_index = 0
        line_ID_index = 0
        load_ID_index = 0
        storage_ID_index = 0

        busbar_sections_IDs = [Bus_obj.BusID]

        # add non-auxilliary buses
        for non_aux_bus_id in 1:N_additional_nonaux_buses
            busID = non_aux_bus_id + grid.N_bus
            switch_id = grid.N_switch+1
            new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(2,busID),Active=true)
            new_bus = Bus(BusID=busID,V_max=Bus_obj.V_max,V_min=Bus_obj.V_min,BusType=0,GeneralSwitch=new_switch)
            push!(grid.Buses,busID => new_bus)
            push!(busbar_sections_IDs,busID)
            grid.N_switch += 1
        end
        grid.N_bus += N_additional_nonaux_buses

        
        # add coupler lines
        switch_id = grid.N_switch+1
        coupler_id = grid.N_branch+1
        status_t = ones(1,grid.N_time_steps)
        new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(2,coupler_id),Active=true,SwitchingStatus_t=status_t)
        new_coupler = Branch(LineID=coupler_id,Fr_bus_ID=Bus_obj.BusID,To_bus_ID=busbar_sections_IDs[2],BranchType=2,GeneralSwitch=new_switch)
        push!(grid.Branches,coupler_id => new_coupler)

        grid.N_switch += 1
        grid.N_branch += 1

        n_gens = length(Bus_obj.ConnectedGensIDs)
        n_lines = length(Bus_obj.ConnectedLinesIDs)
        n_loads = length(Bus_obj.ConnectedLoadsIDs)
        n_strgs = length(Bus_obj.ConnectedStorageIDs)

        gen_bus_interval = 1 : n_gens
        line_bus_interval = n_gens+1 : n_lines+n_gens
        load_bus_interval = n_lines+n_gens+1 : n_lines+n_gens+n_loads
        strg_bus_interval = n_lines+n_gens+n_loads+1 : n_lines+n_gens+n_loads+n_strgs

        # add auxilliary buses
        aux_bus_IDs = []
        for aux_bus_id in 1:N_elements_connected
            switch_id = grid.N_switch+1
            busID = aux_bus_id + grid.N_bus
            new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(2,busID),Active=true)
            new_bus = Bus(BusID=busID,V_max=Bus_obj.V_max,V_min=Bus_obj.V_min,BusType=1,GeneralSwitch=new_switch)
            
            if aux_bus_id in gen_bus_interval

                # add generator auxilliary buses:
                gen_ID_index += 1
                new_bus.ConnectedGensIDs = [Bus_obj.ConnectedGensIDs[gen_ID_index]]
                grid.Generators[Bus_obj.ConnectedGensIDs[gen_ID_index]].GenBus_ID = busID

            elseif aux_bus_id in line_bus_interval
                # add line auxilliary buses
                line_ID_index += 1
                new_bus.ConnectedLinesIDs = [Bus_obj.ConnectedLinesIDs[line_ID_index]]

                if grid.Branches[Bus_obj.ConnectedLinesIDs[line_ID_index]].Fr_bus_ID == Bus_obj.BusID
                    grid.Branches[Bus_obj.ConnectedLinesIDs[line_ID_index]].Fr_bus_ID = busID
                end

                if grid.Branches[Bus_obj.ConnectedLinesIDs[line_ID_index]].To_bus_ID == Bus_obj.BusID
                    grid.Branches[Bus_obj.ConnectedLinesIDs[line_ID_index]].To_bus_ID = busID
                end

            elseif aux_bus_id in load_bus_interval
                # add load auxiliary buses
                load_ID_index += 1
                new_bus.ConnectedLoadsIDs = [Bus_obj.ConnectedLoadsIDs[load_ID_index]]
                grid.Loads[Bus_obj.ConnectedLoadsIDs[load_ID_index]].LoadBus_ID = busID


            elseif aux_bus_id in strg_bus_interval
                # add storage auxiliary buses
                storage_ID_index += 1
                new_bus.ConnectedStorageIDs = [Bus_obj.ConnectedStorageIDs[storage_ID_index]]
                grid.Loads[Bus_obj.ConnectedStorageIDs[storage_ID_index]].StorageBus_ID = busID

            end
            grid.N_switch += 1
            push!(grid.Buses,busID => new_bus)
            push!(aux_bus_IDs,busID)
        end

        grid.N_bus += N_elements_connected

        Bus_obj.ConnectedGensIDs = []
        Bus_obj.ConnectedLinesIDs = []
        Bus_obj.ConnectedLoadsIDs = []
        Bus_obj.ConnectedStorageIDs = []

        # add auxilliary lines
        reconf_aux_lines_IDs = []
        twins = Dict()
        for aux_bus in aux_bus_IDs

            for nonaux_bus in busbar_sections_IDs

                grid.N_switch += 1
                grid.N_branch += 1
                status = nonaux_bus == bus ? 1 : 0
                status_t = ones(1,grid.N_time_steps)*status
                new_switch = Switch(SwitchID=grid.N_switch,Switch_ObjID=(2,grid.N_branch),Active=true,SwitchingStatus=status,SwitchingStatus_t=status_t)
                new_branch = Branch(LineID=grid.N_branch,Fr_bus_ID=nonaux_bus,To_bus_ID=aux_bus,BranchType=1,GeneralSwitch=new_switch)
                push!(grid.Branches,grid.N_branch => new_branch)
                push!(reconf_aux_lines_IDs,grid.N_branch)

                # add the branch ID to both buses
                push!(grid.Buses[aux_bus].ConnectedLinesIDs, grid.N_branch)
                push!(grid.Buses[nonaux_bus].ConnectedLinesIDs, grid.N_branch)

                if nonaux_bus == busbar_sections_IDs[1]
                    push!(twins,grid.N_branch => grid.N_branch+1)
                elseif nonaux_bus == busbar_sections_IDs[2]
                    push!(twins,grid.N_branch => grid.N_branch-1)
                end

                
            end

        end
        # add coupler to the connected lines
        for nonaux_bus in busbar_sections_IDs
            push!(grid.Buses[nonaux_bus].ConnectedLinesIDs, coupler_id)
        end
        
        new_substation = SubStation(SubStationBusID=Bus_obj.BusID,Aux_Buses_IDs=aux_bus_IDs,BusbarSections_IDs=busbar_sections_IDs,
            Reconf_AuxLines_IDs=reconf_aux_lines_IDs,Reconf_CouplerLines_IDs=[coupler_id],Branches_IDs=Bus_obj.ConnectedLinesIDs,
            Generators_IDs=Bus_obj.ConnectedGensIDs,Loads_IDs=Bus_obj.ConnectedLoadsIDs,Storages_IDs=Bus_obj.ConnectedStorageIDs,twins=twins)

        push!(grid.Substations,new_substation.SubStationBusID => new_substation)
        grid.N_substation += 1

        if add_b2b_VSC
            if VSC_capacity isa(Dict)
                add_converter!(grid,busbar_sections_IDs[1],busbar_sections_IDs[2],get(VSC_capacity,bus,1)*grid.S_base;type=:B2B)
            else
                add_converter!(grid,busbar_sections_IDs[1],busbar_sections_IDs[2],VSC_capacity*grid.S_base;type=:B2B)
            end
        end

    end
    Y_bus,b_line = Y_Bus_Grid(grid)
    grid.Y_bus = Y_bus
    grid.b_line = b_line
end

function activate_transmission_switch!(grid ::PowerGrid, branch_ids)
    for id in branch_ids
        grid.Branches[id].GeneralSwitch.Active = true
    end
end

function reduce_grid!(grid ::PowerGrid;clean=false)

    busbar_sections = []
    aux_buses = []
    reconf_lines = []
    couplers = []
    section2substation = Dict()
    for substation_id in keys(grid.Substations)
        non_aux_bus_in_substation = grid.Substations[substation_id].BusbarSections_IDs
        aux_buses_in_substation = grid.Substations[substation_id].Aux_Buses_IDs
        push!(section2substation,non_aux_bus_in_substation[1]=>substation_id)
        push!(section2substation,non_aux_bus_in_substation[2]=>substation_id)
        append!(busbar_sections,non_aux_bus_in_substation)
        append!(aux_buses,aux_buses_in_substation)
        append!(reconf_lines,grid.Substations[substation_id].Reconf_AuxLines_IDs)
        append!(couplers,grid.Substations[substation_id].Reconf_CouplerLines_IDs)
    end

    # merge connected buses inside substations
    for substation_id in keys(grid.Substations)
        reconf_line_ids = grid.Substations[substation_id].Reconf_AuxLines_IDs
        coupler_line_ids = grid.Substations[substation_id].Reconf_CouplerLines_IDs
        removed_lines = []
        # check reconfiguration lines
        for reconf_id in reconf_line_ids
            if reconf_id ∉ removed_lines
                reconf_line = grid.Branches[reconf_id]
                #check if reconf line is on to merge buses
                if reconf_line.GeneralSwitch.SwitchingStatus == 1 
                    Fr_bus_id = reconf_line.Fr_bus_ID # nonaux root bus (will not be removed)
                    To_bus_id = reconf_line.To_bus_ID # aux bus (will be merged)
                
                    line_ids = []
                    for line_id in grid.Buses[To_bus_id].ConnectedLinesIDs
                        if grid.Branches[line_id].BranchType == 0
                            append!(line_ids,line_id)
                        end
                    end
                    load_ids = grid.Buses[To_bus_id].ConnectedLoadsIDs
                    gen_ids = grid.Buses[To_bus_id].ConnectedGensIDs
                    storage_ids = grid.Buses[To_bus_id].ConnectedStorageIDs

                    append!(grid.Buses[Fr_bus_id].ConnectedLinesIDs,line_ids)
                    append!(grid.Buses[Fr_bus_id].ConnectedLoadsIDs,load_ids)
                    append!(grid.Buses[Fr_bus_id].ConnectedGensIDs,gen_ids)
                    append!(grid.Buses[Fr_bus_id].ConnectedStorageIDs,storage_ids)

                    [grid.Generators[gen_id].GenBus_ID = Fr_bus_id for gen_id in gen_ids]
                    [grid.Loads[load_id].LoadBus_ID = Fr_bus_id for load_id in load_ids]
                    [grid.Storages[storage_id].StorageBus_ID = Fr_bus_id for storage_id in storage_ids]
                    
                    for branch_id in line_ids
                        branch_obj = grid.Branches[branch_id]
                        if branch_obj.Fr_bus_ID == To_bus_id
                            grid.Branches[branch_id].Fr_bus_ID = Fr_bus_id
                        else
                            grid.Branches[branch_id].To_bus_ID = Fr_bus_id
                        end
                    end
                    push!(removed_lines,reconf_id)
                    index_at_nonaux = findfirst(x->x==reconf_id,grid.Buses[Fr_bus_id].ConnectedLinesIDs)
                    deleteat!(grid.Buses[Fr_bus_id].ConnectedLinesIDs,index_at_nonaux)
                    delete!(grid.Branches,reconf_id)
                    delete!(grid.Buses,To_bus_id)
                    grid.N_bus -= 1
                    grid.N_branch -= 1
                    grid.N_switch -= 2
                else
                    Fr_bus = reconf_line.Fr_bus_ID # nonaux root bus (will not be removed)
                    To_bus = reconf_line.To_bus_ID # aux bus (will be merged)

                    index_at_Fr = findfirst(x-> x==reconf_id, grid.Buses[Fr_bus].ConnectedLinesIDs)
                    deleteat!(grid.Buses[Fr_bus].ConnectedLinesIDs, index_at_Fr)
                    if To_bus in keys(grid.Buses)
                        index_at_To = findfirst(x-> x==reconf_id, grid.Buses[To_bus].ConnectedLinesIDs)
                        deleteat!(grid.Buses[To_bus].ConnectedLinesIDs, index_at_To)
                    end
                end

            end
        end
        if clean
            # check couplers
            for coupler_id in coupler_line_ids
                coupler_branch = grid.Branches[coupler_id]
                if coupler_branch.GeneralSwitch.SwitchingStatus == 1
                    transferable_branches = collect(setdiff(Set(grid.Buses[coupler_branch.To_bus_ID].ConnectedLinesIDs),Set(removed_lines)))
                    append!(grid.Buses[coupler_branch.Fr_bus_ID].ConnectedGensIDs,grid.Buses[coupler_branch.To_bus_ID].ConnectedGensIDs)
                    append!(grid.Buses[coupler_branch.Fr_bus_ID].ConnectedLinesIDs,transferable_branches)
                    append!(grid.Buses[coupler_branch.Fr_bus_ID].ConnectedLoadsIDs,grid.Buses[coupler_branch.To_bus_ID].ConnectedLoadsIDs)
                    append!(grid.Buses[coupler_branch.Fr_bus_ID].ConnectedStorageIDs,grid.Buses[coupler_branch.To_bus_ID].ConnectedStorageIDs)
                    
                    index_at_bus = findfirst(x->x==coupler_id, grid.Buses[coupler_branch.Fr_bus_ID].ConnectedLinesIDs)
                    if index_at_bus !== nothing
                        deleteat!(grid.Buses[coupler_branch.Fr_bus_ID].ConnectedLinesIDs,index_at_bus)
                    end

                    [grid.Generators[gen_id].GenBus_ID = coupler_branch.Fr_bus_ID for gen_id in grid.Buses[coupler_branch.To_bus_ID].ConnectedGensIDs]
                    [grid.Loads[load_id].LoadBus_ID = coupler_branch.Fr_bus_ID for load_id in grid.Buses[coupler_branch.To_bus_ID].ConnectedLoadsIDs]
                    [grid.Storages[storage_id].StorageBus_ID = coupler_branch.Fr_bus_ID for storage_id in grid.Buses[coupler_branch.To_bus_ID].ConnectedStorageIDs]
                    
                    for branch_id in transferable_branches
                        branch_obj = grid.Branches[branch_id]
                        if branch_obj.Fr_bus_ID == coupler_branch.To_bus_ID
                            grid.Branches[branch_id].Fr_bus_ID = coupler_branch.Fr_bus_ID
                        else
                            grid.Branches[branch_id].To_bus_ID = coupler_branch.Fr_bus_ID
                        end
                    end
                    grid.Buses[coupler_branch.To_bus_ID].GeneralSwitch.SwitchingStatus = 0 
                end
                

            end
        end
    end

    if clean
        for branch in keys(grid.Branches)
            if grid.Branches[branch].BranchType != 0
                delete!(grid.Branches,branch)
                grid.N_branch -= 1
                grid.N_switch -= 1
            end
        end
    else
        for branch in keys(grid.Branches)
            if grid.Branches[branch].BranchType == 1
                delete!(grid.Branches,branch)
                grid.N_branch -= 1
                grid.N_switch -= 1
            end
        end
    end

    # remove inactive buses
    for bus_id in keys(grid.Buses)
        if grid.Buses[bus_id].GeneralSwitch.SwitchingStatus == 0
            delete!(grid.Buses,bus_id)
            grid.N_bus -= 1
            grid.N_switch -= 1
        end
    end

    if clean
        for bus_id in keys(grid.Buses)
            if bus_id > grid.N_bus
                _rename_bus!(grid, bus_id)
            end
        end
    else
        # TODO: probably something should be done here
        for bus_id in keys(grid.Buses)
            if bus_id > grid.N_bus
                _rename_bus!(grid, bus_id)
            end
        end
    end
    
end

function _rename_bus!(grid ::PowerGrid,bus_id)
    bus_substation = find_substation(grid,bus_id)
    if bus_substation ∉ keys(grid.Buses)
        bus_obj = grid.Buses[bus_id]
        gen_ids =  bus_obj.ConnectedGensIDs
        line_ids = bus_obj.ConnectedLinesIDs
        load_ids = bus_obj.ConnectedLoadsIDs
        storage_ids = bus_obj.ConnectedStorageIDs
        new_bus = deepcopy(grid.Buses[bus_id])
        new_bus.BusID = bus_substation

        [grid.Generators[g].GenBus_ID = bus_substation for g in gen_ids]
        [grid.Loads[d].LoadBus_ID = bus_substation for d in load_ids]
        [grid.Storages[s].StorageBus_ID = bus_substation for s in storage_ids]

        for line in line_ids
            if grid.Branches[line].Fr_bus_ID == bus_id
                grid.Branches[line].Fr_bus_ID = bus_substation
            else
                grid.Branches[line].To_bus_ID = bus_substation
            end
        end
        delete!(grid.Buses,bus_id)
        push!(grid.Buses, bus_substation => new_bus)

        index = findfirst(x->x==bus_id,grid.Substations[bus_substation].BusbarSections_IDs)
        grid.Substations[bus_substation].BusbarSections_IDs[index] = bus_substation
    else
        # do something
        new_bus_id = maximum([key for key in keys(grid.Buses) if key ≤ grid.N_bus])+1
        bus_obj = grid.Buses[bus_id]
        gen_ids =  bus_obj.ConnectedGensIDs
        line_ids = bus_obj.ConnectedLinesIDs
        load_ids = bus_obj.ConnectedLoadsIDs
        storage_ids = bus_obj.ConnectedStorageIDs
        new_bus = deepcopy(grid.Buses[bus_id])
        new_bus.BusID = new_bus_id

        [grid.Generators[g].GenBus_ID = new_bus_id for g in gen_ids]
        [grid.Loads[d].LoadBus_ID = new_bus_id for d in load_ids]
        [grid.Storages[s].StorageBus_ID = new_bus_id for s in storage_ids]

        for line in line_ids
            if grid.Branches[line].Fr_bus_ID == bus_id
                grid.Branches[line].Fr_bus_ID = new_bus_id
            else
                grid.Branches[line].To_bus_ID = new_bus_id
            end
        end
        grid.Substations[bus_substation].BusbarSections_IDs[2] = new_bus_id
        delete!(grid.Buses,bus_id)
        push!(grid.Buses, new_bus_id => new_bus)
    end
end

function find_substation(grid ::PowerGrid,bus_id)
    return [s for s in keys(grid.Substations) if bus_id in grid.Substations[s].BusbarSections_IDs][1]
end

function _get_edge_list(grid ::PowerGrid,from_line_loading=true)
    edge_list = []

    if from_line_loading
        for row in eachrow(grid.LineLoading)
            Fr_bus = row[:FromBus]
            To_bus = row[:ToBus]
            push!(edge_list,(Fr_bus,To_bus))
        end

    else
        for branch_id in keys(grid.Branches)
            Fr_bus = grid.Branches[branch_id].Fr_bus_ID
            To_bus = grid.Branches[branch_id].To_bus_ID
            push!(edge_list,(Fr_bus,To_bus))
        end
    end

    return edge_list
end

function apply_switching_status!(grid ::PowerGrid)

    for line_id in keys(grid.Z_lines)
        grid.Branches[line_id].GeneralSwitch.SwitchingStatus = grid.Z_lines[line_id]
    end
    
end

