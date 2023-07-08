function parse_matpower_case(m_file_path; start_node_from_1=true)
    bus_data_raw = [] #	bus_i	type	Pd	Qd	Gs	Bs	area	Vm      Va	baseKV	zone	Vmax	Vmin
    gen_data_raw = [] #	bus	Pg      Qg	Qmax	Qmin	Vg	mBase       status	Pmax	Pmin	pc1 pc2 qlcmin qlcmax qc2min qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
    branch_data_raw = [] #	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle   status angmin angmax
    dc_bus_data_raw = [] #  busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
    conv_data_raw = [] #  busdc_i busac_i type_dc type_ac P_g   Q_g islcc  Vtar    rtf xtf  transformer tm   bf filter    rc      xc  reactor   basekVac    Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop      Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin
    dc_branch_data_raw = [] #  fbusdc  tbusdc  r      l        c   rateA   rateB   rateC   status
    gen_cost_data_raw = [] #   2	startup	shutdown	n	c(n-1)	...	c0
	dc_link_data_raw = [] # fbus	tbus	status	Pf	Pt	Qf	Qt	Vf	Vt	Pmin	Pmax	QminF	QmaxF	QminT	QmaxT	loss0	loss1

    open(m_file_path) do f
        read_bus_data = false
        read_gen_data = false
        read_branch_data = false
        read_dc_bus_data = false
        read_conv_data = false
        read_dc_branch_data = false
        read_gen_cost_data = false
        read_dc_link_data = false
        
        while !eof(f)
            x = readline(f)

            
            if startswith(x, "%")
                    continue
            end

            if occursin("mpc.bus ",x)
                    read_bus_data = true
                    continue
            end
            if occursin("mpc.gen ",x)
                    read_gen_data = true
                    continue
            end
            if occursin("mpc.branch ",x)
                    read_branch_data = true
                    continue
            end
            if occursin("mpc.busdc",x)
                    read_dc_bus_data = true
                    continue
            end
            if occursin("mpc.convdc",x)
                    read_conv_data = true
                    continue
            end
            if occursin("mpc.branchdc",x)
                    read_dc_branch_data = true
                    continue
            end
            if occursin("mpc.dcline ",x)
                    read_dc_link_data = true
                    continue
            end
            if occursin("mpc.gencost",x)
                    read_gen_cost_data = true
                    continue
            end

            if occursin("]",x) # close flages
                    read_bus_data = false
                    read_gen_data = false
                    read_branch_data = false
                    read_dc_bus_data = false
                    read_conv_data = false
                    read_dc_branch_data = false
                    read_gen_cost_data = false
		    read_dc_link_data = false
                    continue
            end

            if read_bus_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(bus_data_raw,b)
                    continue
            end

            if read_gen_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(gen_data_raw,b)
                    continue
            end

            if read_branch_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(branch_data_raw,b)
                    continue
            end

            if read_dc_bus_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(dc_bus_data_raw,b)
                    continue
            end

            if read_conv_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(conv_data_raw,b)
                    continue
            end

            if read_dc_branch_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(dc_branch_data_raw,b)
                    continue
            end

            if read_dc_link_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(dc_link_data_raw,b)
                    continue
            end

            if read_gen_cost_data
                    parts = split(x, "%", limit=2)
                    x = strip(parts[1])
                    a = split(x)
                    a[end] = replace(a[end],";"=>"")
                    b = map(γ->parse(Float64,γ),a)
                    push!(gen_cost_data_raw,b)
                    continue
            end

        end
    end

    return _create_grid_from_raw_data(bus_data_raw,gen_data_raw,branch_data_raw,dc_bus_data_raw,conv_data_raw,dc_branch_data_raw,gen_cost_data_raw,dc_link_data_raw;start_node_from_1=start_node_from_1)

end

function _create_grid_from_raw_data(bus_data_raw,gen_data_raw,branch_data_raw,dc_bus_data_raw,conv_data_raw,dc_branch_data_raw,gen_cost_data_raw,dc_link_data_raw;start_node_from_1=true)

    grid = PowerGrid()
    if start_node_from_1
        ac_node_map = Dict()
        dc_node_map = Dict()
        ac_node_id_in_grid = 0
        dc_node_id_in_grid = 0
    end
    # add Buses and loads
    for row in eachrow(bus_data_raw)
        data = row[1] #	bus_i	type	Pd	Qd	Gs	Bs	area	Vm      Va	baseKV	zone	Vmax	Vmin
        if start_node_from_1
                ac_node_id_in_grid += 1
                push!(ac_node_map, Int64(data[1]) => ac_node_id_in_grid)
                add_bus!(grid,ac_node_id_in_grid, V_max=data[12],V_min=data[13],δ_max=0.6,δ_min=-0.6)
                add_load!(grid,ac_node_id_in_grid, data[3],data[4])
        else
                add_bus!(grid,Int64(data[1]), V_max=data[12],V_min=data[13],δ_max=0.6,δ_min=-0.6)
                add_load!(grid,Int64(data[1]), data[3],data[4])
        end 
        # grid.Loads[grid.N_load].Pd_t = [data[3]]
    end
    # add DC Buses
    for row in eachrow(dc_bus_data_raw)
        data = row[1] #  busdc_i grid    Pdc     Vdc     basekVdc    Vdcmax  Vdcmin  Cdc
        if start_node_from_1
                dc_node_id_in_grid += 1
                push!(dc_node_map, Int64(data[1]) => dc_node_id_in_grid)
                add_dc_bus!(grid,dc_node_id_in_grid,V_max=data[6],V_min=data[7])
        else
                add_dc_bus!(grid,Int64(data[1]),V_max=data[6],V_min=data[7])
        end
        
    end

    # add branches 
    for row in eachrow(branch_data_raw)
        data = row[1] #	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle   status angmin angmax
        if start_node_from_1
                fbus = ac_node_map[Int64(data[1])]
                tbus = ac_node_map[Int64(data[2])]
                add_branch!(grid,fbus,tbus;r_ohms=data[3],x_ohms=data[4],b_ohms=data[5],rating_pu=data[6]/100,branch_type=0)
        else
                add_branch!(grid,Int64(data[1]),Int64(data[2]);r_ohms=data[3],x_ohms=data[4],b_ohms=data[5],rating_pu=data[6]/100,branch_type=0)
        end
    end

    #add DC Branches
    for row in  eachrow(dc_branch_data_raw)
        data = row[1] #  fbusdc  tbusdc  r      l        c   rateA   rateB   rateC   status
        if start_node_from_1
                fbus = dc_node_map[Int64(data[1])]
                tbus = dc_node_map[Int64(data[2])]
                add_dc_branch!(grid,fbus,tbus;r_ohms=data[3],x_ohms=data[4],b_ohms=data[5],rating_pu=data[6]/100)
        else
                add_dc_branch!(grid,Int64(data[1]),Int64(data[2]);r_ohms=data[3],x_ohms=data[4],b_ohms=data[5],rating_pu=data[6]/100)
        end
    end

    # add Generators
    for row in zip(gen_data_raw,gen_cost_data_raw)
        data = row[1] #	bus	Pg      Qg	Qmax	Qmin	Vg	mBase       status	Pmax	Pmin	pc1 pc2 qlcmin qlcmax qc2min qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
        cost = row[2] #   2	startup	shutdown	n	c(n-1)	...	c0
        if start_node_from_1
                g_bus = ac_node_map[Int64(data[1])]
                add_generator!(grid,g_bus,cost[end],cost[end-1],cost[end-2],data[9],data[10],data[4],data[5]; GenType=:base_type,
                        start_up_cost=cost[2],shut_down_cost=cost[3],min_up_time=1,min_down_time=1,Δ_up=0,Δ_down=0,Pg=0,Qg=0)
        else
                add_generator!(grid,Int64(data[1]),cost[end],cost[end-1],cost[end-2],data[9],data[10],data[4],data[5]; GenType=:base_type,
                        start_up_cost=cost[2],shut_down_cost=cost[3],min_up_time=1,min_down_time=1,Δ_up=0,Δ_down=0,Pg=0,Qg=0)
        end
    end

    # add Converters
    for row in eachrow(conv_data_raw)
        data = row[1] #  busdc_i busac_i type_dc type_ac P_g   Q_g islcc  Vtar    rtf xtf  transformer tm   bf filter    rc      xc  reactor   basekVac    Vmmax   Vmmin   Imax    status   LossA LossB  LossCrec LossCinv  droop      Pdcset    Vdcset  dVdcset Pacmax Pacmin Qacmax Qacmin
        if start_node_from_1
                dc_bus = dc_node_map[Int64(data[1])]
                ac_bus = ac_node_map[Int64(data[2])]
                add_converter!(grid,dc_bus,ac_bus,data[end-3];type=:ACDC,loss_a=data[23], loss_b=data[24])
        else
                add_converter!(grid,Int64(data[1]),Int64(data[2]),data[end-3];type=:ACDC,loss_a=data[23], loss_b=data[24])
        end
    end

     # add DC links
     for row in eachrow(dc_link_data_raw)
        data = row[1] # fbus	tbus	status	Pf	Pt	Qf	Qt	Vf	Vt	Pmin	Pmax	QminF	QmaxF	QminT	QmaxT	loss0	loss1
        if start_node_from_1
                fbus = ac_node_map[Int64(data[1])]
                tbus = ac_node_map[Int64(data[2])]
                add_dc_link!(grid,fbus,tbus,data[11])
        else
                add_dc_link!(grid,Int64(data[1]),Int64(data[2]),data[11])
        end
     end

    Y_bus,b_line = Y_Bus_Grid(grid)
    grid.Y_bus = Y_bus
    grid.b_line = b_line

    return grid
end
