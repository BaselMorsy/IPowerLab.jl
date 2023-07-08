# include("../Components/Grid.jl")

using JLD2
using FileIO

function add_bus!(grid ::PowerGrid,bus_id;V_max=1.1,V_min=0.9,δ_max=0.6,δ_min=-0.6)
    switch_id = grid.N_switch + 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(1,bus_id))
    new_bus = Bus(BusID=bus_id,GeneralSwitch=new_switch,V_max=V_max,V_min=V_min,δ_max=δ_max,δ_min=δ_min)
    push!(grid.Buses,bus_id => new_bus)
    grid.N_bus += 1
    grid.N_switch += 1
end

function add_dc_bus!(grid ::PowerGrid,bus_id;V_max=1.1,V_min=0.9)
    switch_id = grid.N_switch + 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(9,bus_id))
    new_bus = DCBus(BusID=bus_id,GeneralSwitch=new_switch,V_max=V_max,V_min=V_min)
    push!(grid.DCBuses,bus_id => new_bus)
    grid.N_dc_bus += 1
    grid.N_switch += 1
end

function add_generator!(grid ::PowerGrid,gen_bus,C0,C1,C2,Pg_max,Pg_min,Qg_max,Qg_min; GenType=:base_type,
        start_up_cost=0,shut_down_cost=0,min_up_time=1,min_down_time=1,Δ_up=0,Δ_down=0,Pg=0,Qg=0,gen_id=nothing,bus_type=:AC)

    if gen_id === nothing
        if bus_type == :AC
            gen_id = grid.N_gen + 1
            grid.N_gen += 1
        elseif bus_type == :DC
            gen_id = grid.N_dc_gen + 1
            grid.N_dc_gen += 1
        else
            error("Unrecognized `bus_type`: ($bus_type)")
        end
    else
        if bus_type == :AC
            grid.N_gen += 1
        elseif bus_type == :DC
            grid.N_dc_gen += 1
        else
            error("Unrecognized `bus_type`: ($bus_type)")
        end
    end

    switch_id = grid.N_switch + 1
    grid.N_switch += 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(3,gen_id))
    new_gen = Generator(GenID=gen_id,GenBus_ID=gen_bus,GeneralSwitch=new_switch,C0=C0,C1=C1,C2=C2,
        Pg_max=Pg_max,Pg_min=Pg_min,Qg_max=Qg_max,Qg_min=Qg_min,GenType=GenType,start_up_cost=start_up_cost,
        shut_down_cost=shut_down_cost,min_up_time=min_up_time,min_down_time=min_down_time,Δ_up=Δ_up,Δ_down=Δ_down,
        Pg=Pg,Qg=Qg)

    if bus_type == :AC
        push!(grid.Buses[gen_bus].ConnectedGensIDs,gen_id)
        push!(grid.Generators, gen_id => new_gen)
    elseif bus_type == :DC
        push!(grid.DCBuses[gen_bus].ConnectedGensIDs,gen_id)
        push!(grid.DCGenerators, gen_id => new_gen)
    end
    
end

function add_load!(grid ::PowerGrid,load_bus,Pd,Qd;Pd_profile_scaled=[],Qd_profile_scaled=[],Shedding_Cost=1000,load_id=nothing,bus_type=:AC)

    if load_id === nothing
        if bus_type == :AC
            load_id = grid.N_load + 1
            grid.N_load += 1
        elseif bus_type == :DC
            load_id = grid.N_dc_load + 1
            grid.N_dc_load += 1
        else
            error("Unrecognized `bus_type`: ($bus_type)")
        end
    else
        if bus_type == :AC
            grid.N_load += 1
        elseif bus_type == :DC
            grid.N_dc_load += 1
        else
            error("Unrecognized `bus_type`: ($bus_type)")
        end
    end
    switch_id = grid.N_switch + 1
    grid.N_switch += 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(4,load_id))
    new_load = Load(LoadID=load_id,LoadBus_ID=load_bus,Pd=Pd,Qd=Qd,Pd_profile_scaled=Pd_profile_scaled,Qd_profile_scaled=Qd_profile_scaled,
        Shedding_Cost=Shedding_Cost,GeneralSwitch=new_switch,Pd_t=[Pd], Qd_t=[Qd])
    if bus_type == :AC
        push!(grid.Loads, load_id => new_load)
        push!(grid.Buses[load_bus].ConnectedLoadsIDs, load_id)
    elseif bus_type == :DC
        push!(grid.DCLoads, load_id => new_load)
        push!(grid.DCBuses[load_bus].ConnectedLoadsIDs, load_id)
    end
    
end

function add_branch!(grid ::PowerGrid,Fr_bus_ID,To_bus_ID;r_ohms=Inf,x_ohms=Inf,b_ohms=0,rating_pu=50,branch_type=0,branch_id=nothing)

    if branch_id === nothing
        branch_id = grid.N_branch + 1
    end
    switch_id = grid.N_switch + 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(2,branch_id))
    new_branch = Branch(LineID=branch_id,Fr_bus_ID=Fr_bus_ID,To_bus_ID=To_bus_ID,GeneralSwitch=new_switch,
        BranchType=branch_type,r=r_ohms,x=x_ohms,b=b_ohms,rating=rating_pu)
    push!(grid.Branches, branch_id => new_branch)
    grid.N_branch += 1
    grid.N_switch += 1
    push!(grid.Buses[Fr_bus_ID].ConnectedLinesIDs, branch_id)
    push!(grid.Buses[To_bus_ID].ConnectedLinesIDs, branch_id)


end

function add_dc_branch!(grid ::PowerGrid,Fr_bus_ID,To_bus_ID;r_ohms=Inf,x_ohms=Inf,b_ohms=0,rating_pu=50,branch_id=nothing)

    if branch_id === nothing
        branch_id = grid.N_dc_branch + 1
    end
    switch_id = grid.N_switch + 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(10,branch_id))
    new_branch = Branch(LineID=branch_id,Fr_bus_ID=Fr_bus_ID,To_bus_ID=To_bus_ID,GeneralSwitch=new_switch,
        r=r_ohms,x=x_ohms,b=b_ohms,rating=rating_pu)
    push!(grid.DCBranches, branch_id => new_branch)
    grid.N_dc_branch += 1
    grid.N_switch += 1
    push!(grid.DCBuses[Fr_bus_ID].ConnectedLinesIDs, branch_id)
    push!(grid.DCBuses[To_bus_ID].ConnectedLinesIDs, branch_id)


end

function add_converter!(grid ::PowerGrid,DC_Bus_ID,AC_Bus_ID,rate;type=:ACDC,loss_a=0,loss_b=0,conv_id=nothing)
    if conv_id === nothing
        conv_id = grid.N_converter + 1
    end
    grid.N_converter += 1
    switch_id = grid.N_switch + 1
    grid.N_switch += 1
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(8,conv_id))
    new_converter = Converter(Conv_ID=conv_id,DC_Bus_ID=DC_Bus_ID,AC_Bus_ID=AC_Bus_ID,rate=rate,GeneralSwitch=new_switch,type=type,loss_a=loss_a,loss_b=loss_b)

    if type == :ACDC
        add_generator!(grid,DC_Bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:DC)
        new_converter.gen_dc_id = grid.N_dc_gen
        add_generator!(grid,AC_Bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:AC)
        new_converter.gen_ac_id = grid.N_gen
    elseif type == :B2B
        add_generator!(grid,DC_Bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:AC)
        new_converter.gen_dc_id = grid.N_gen
        add_generator!(grid,AC_Bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:AC)
        new_converter.gen_ac_id = grid.N_gen
    end
    push!(grid.Converters, conv_id => new_converter)
end

function add_dc_link!(grid ::PowerGrid,Fr_bus_ID,To_bus_ID,rate)
    grid.N_dc_links += 1
    grid.N_switch += 1
    dc_link_id = grid.N_dc_links
    switch_id = grid.N_switch
    new_switch = Switch(SwitchID=switch_id,Switch_ObjID=(11,dc_link_id))
    new_link = DCLink(LineID=dc_link_id,Fr_bus_ID=Fr_bus_ID,To_bus_ID=To_bus_ID,rating=rate,GeneralSwitch=new_switch)
    add_generator!(grid,Fr_bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:AC)
    new_link.Fr_gen_ID = grid.N_gen
    add_generator!(grid,To_bus_ID,0,0,0,rate,-rate,rate,-rate,GenType=:virtual,bus_type=:AC)
    new_link.To_gen_ID = grid.N_gen
    push!(grid.DCLinks, dc_link_id => new_link)
end

function split_element!(grid ::PowerGrid, object_type ::Symbol, object_id ::Int64;n=2,modularization=:discrete)
    
    if object_type == :Converter
        my_converter = grid.Converters[object_id]
        DC_Bus_ID = my_converter.DC_Bus_ID
        AC_Bus_ID = my_converter.AC_Bus_ID
        rate = my_converter.rate
        type = my_converter.type
        if modularization == :discrete
            grid.Converters[object_id].rate *= 1/n
            grid.Generators[grid.Converters[object_id].gen_dc_id].Pg_max *= 1/n
            grid.Generators[grid.Converters[object_id].gen_ac_id].Pg_max *= 1/n
            grid.Generators[grid.Converters[object_id].gen_dc_id].Pg_min *= 1/n
            grid.Generators[grid.Converters[object_id].gen_ac_id].Pg_min *= 1/n
            for i in collect(1:n)
                add_converter!(grid,DC_Bus_ID,AC_Bus_ID,rate*1/n,type=type)
            end
        elseif modularization == :continuous
            if n > 2
                error("Continous capacity applies only for 2 converters!")
            else
                add_converter!(grid,DC_Bus_ID,AC_Bus_ID,rate,type=type)
                grid.N_conv_duplets += 1
                push!(grid.Converter_Duplets, grid.N_conv_duplets => (my_converter.Conv_ID, grid.N_converter))
            end

        else
            error("Unknown converter modularization type ($modularization)")
        end
        
    elseif object_type == :Generator
        my_gen = deepcopy(grid.Generators[object_id])
        add_generator!(grid, my_gen.GenBus_ID, my_gen.C0, my_gen.C1, my_gen.C2, my_gen.Pg_max, my_gen.Pg_min, my_gen.Qg_max, my_gen.Qg_min;
            Δ_up=my_gen.Δ_up, Δ_down=my_gen.Δ_down)
        grid.N_gen_duplets += 1
        push!(grid.Generator_Duplets, grid.N_gen_duplets => (my_gen.GenID, grid.N_gen))      
    elseif object_type == :Branch
    elseif object_type == :DCBranch
    elseif object_type == :DCLink
    elseif object_type == :Load
    else
        error("Object type ($object_type) is not splittable.")
    end
end

function merge_element!(grid ::PowerGrid, object_type ::Symbol, duplete_id ::Int64)
    if object_type == :Converter

        root_id = grid.Converter_Duplets[duplete_id][1]
        duplicate_id = grid.Converter_Duplets[duplete_id][2]

        root_gen_ac_id = grid.Converters[root_id].gen_ac_id
        root_gen_dc_id = grid.Converters[root_id].gen_dc_id

        duplicate_gen_ac_id = grid.Converters[duplicate_id].gen_ac_id
        duplicate_gen_dc_id = grid.Converters[duplicate_id].gen_dc_id

        for t in keys(grid.Generators[root_gen_ac_id].Pg_tk)
            for k in keys(grid.Generators[root_gen_ac_id].Pg_tk[t])
                grid.Generators[root_gen_ac_id].Pg_tk[t][k] += deepcopy(grid.Generators[duplicate_gen_ac_id].Pg_tk[t][k])
                grid.DCGenerators[root_gen_dc_id].Pg_tk[t][k] += deepcopy(grid.DCGenerators[duplicate_gen_dc_id].Pg_tk[t][k])
            end
        end

        dublicate_ac_gen_bus = grid.Generators[duplicate_gen_ac_id].GenBus_ID
        dublicate_dc_gen_bus = grid.Generators[duplicate_gen_dc_id].GenBus_ID

        deleteat!(grid.Buses[dublicate_ac_gen_bus].ConnectedGensIDs, findall(x->x==duplicate_gen_ac_id, grid.Buses[dublicate_ac_gen_bus].ConnectedGensIDs))
        deleteat!(grid.DCBuses[dublicate_dc_gen_bus].ConnectedGensIDs, findall(x->x==duplicate_gen_dc_id, grid.Buses[dublicate_dc_gen_bus].ConnectedGensIDs))

        delete!(grid.Generators, duplicate_gen_ac_id)
        delete!(grid.DGenerators, duplicate_gen_dc_id)
        delete!(grid.Converters, duplicate_id)
        grid.N_gen -= 1
        grid.N_dc_gen -= 1
        grid.N_converter -= 1
        grid.N_switch -= 3
        delete!(grid.Converter_Duplets, duplete_id)
        
    elseif object_type == :Generator

        root_id = grid.Generator_Duplets[duplete_id][1]
        duplicate_id = grid.Generator_Duplets[duplete_id][2]

        for t in keys(grid.Generators[root_id].Pg_tk)
            for k in keys(grid.Generators[root_id].Pg_tk[t])
                grid.Generators[root_id].Pg_tk[t][k] += deepcopy(grid.Generators[duplicate_id].Pg_tk[t][k])
            end
        end

        dublicate_gen_bus = grid.Generators[duplicate_id].GenBus_ID
        deleteat!(grid.Buses[dublicate_gen_bus].ConnectedGensIDs, findall(x->x==duplicate_id, grid.Buses[dublicate_gen_bus].ConnectedGensIDs))
        delete!(grid.Generators, duplicate_id)
        grid.N_gen -= 1
        grid.N_switch -= 1
        grid.N_gen_duplets -= 1
        delete!(grid.Generator_Duplets, duplete_id)

    elseif object_type == :Branch
    elseif object_type == :DCBranch
    elseif object_type == :DCLink
    elseif object_type == :Load
    else
        error("Object type ($object_type) is not splittable.")
    end
end

function add_storage!(grid ::PowerGrid)
    # To be implemented    
end

function apply_single_load_profile!(grid ::PowerGrid, load_profile)
    for load in keys(grid.Loads)
        grid.Loads[load].Pd_t = grid.Loads[load].Pd*load_profile
        grid.Loads[load].Qd_t = grid.Loads[load].Qd*load_profile
    end
end

function load_custom_grid(CASE_DIR;S_base=100)
    bus_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","bus.csv"),DataFrame ;copycols = true))
    line_data_raw = DataFrame(CSV.read(string(CASE_DIR,"/","branch.csv"),DataFrame ;copycols = true))
    gen_data_raw =  DataFrame(CSV.read(string(CASE_DIR,"/","gen.csv"),DataFrame ;copycols = true))

    grid = PowerGrid(GridID=1)
    grid.S_base = S_base

    for bus in eachrow(bus_data_raw)
        add_bus!(grid,bus[:bus_i],V_max=bus[:Vmax],V_min=bus[:Vmin],δ_max=bus[:delta_max],δ_min=bus[:delta_min])
    end

    for load in eachrow(bus_data_raw)
        add_load!(grid,load[:bus_i],load[:Pd],load[:Qd])
    end

    for gen in eachrow(gen_data_raw)
        add_generator!(grid,gen[:bus],gen[:C0],gen[:C1],gen[:C2],gen[:Pmax],gen[:Pmin],gen[:Qmax],gen[:Qmin],
            start_up_cost=gen[:C_SU],shut_down_cost=gen[:C_SD],min_up_time=gen[:min_up_time],min_down_time=gen[:min_down_time],Δ_up=gen[:ramp_up],
            Δ_down=gen[:ramp_down],Pg=gen[:Pg],Qg=gen[:Qg])
    end

    for line in eachrow(line_data_raw)
        add_branch!(grid,line[:fbus],line[:tbus],r_ohms=line[:r],x_ohms=line[:x],b_ohms=line[:b],rating_pu=line[:rating]./S_base)
    end

    Y_bus,b_line = Y_Bus_Grid(grid)
    grid.Y_bus = Y_bus
    grid.b_line = b_line

    return grid
    
end

function load_grid_snapshot(grid ::PowerGrid,t)
    # grid2load = deepcopy(grid)
end

function save_grid_object(grid ::PowerGrid,save_path)
    save(save_path*".jld2","mygrid",grid)
end

function load_grid_object(load_path)
    return load(load_path*".jld2")["mygrid"]
end

function idealize_branch!(grid ::PowerGrid, branch_id)
    branch_prop = grid.Branches[branch_id]
    value = pop!(grid.Branches, branch_id, nothing) 
       
    add_dc_link!(grid, branch_prop.Fr_bus_ID, branch_prop.To_bus_ID, branch_prop.rating)
    link_id = grid.N_dc_links
    push!(grid.Idealized_Branches_Mapping, branch_prop.LineID => link_id)
    push!(grid.Idealized_Branches, branch_id => branch_prop)

end

function idealize_substation(grid::PowerGrid, substation_ids)
    relaxed_nodes = []
    relaxed_lines = []
    for substation_id in substation_ids
        connected_lines = grid.Buses[substation_id].ConnectedLinesIDs
        append!(relaxed_nodes, substation_id)
        for line_id in connected_lines
            append!(relaxed_lines, line_id)
        end
    end
    return relaxed_nodes, relaxed_lines
end

function deidealize_substation!(grid ::PowerGrid, substation_ids)

end