function _get_links(url::AbstractString)
    response = HTTP.get(url)
    html = String(response.body)

    regex = Regex("<a\\s+(?:[^>]*?\\s+)?href=\"([^\"]*)\"")
    matches = eachmatch(regex, html)

    links = []
    for match in matches
        push!(links, match.captures[1])
    end

    return links
end

function download_available_UC_cases_on_date(date="2017-02-18")
    url = "https://axavier.org/UnitCommitment.jl/0.3/instances/matpower"  # Replace with the actual URL
    links = _get_links(url)
    case_names = links[11:end-1]
    for case in case_names
        download_UC_case(case)
    end
end

function download_UC_case(case_name::AbstractString; date="2017-02-18")
    """
    case_name: case name from "https://axavier.org/UnitCommitment.jl/0.3/instances/matpower" -> has to be a matpower case if you're 
                interested in power flows
    date: anything between "2017-01-01" and "2017-12-31", default is 18th of Feb. because it's my birthday :D 
    """
    INSTANCES_URL = "https://axavier.org/UnitCommitment.jl/0.3/instances"
    name = "/matpower/$case_name/$date"
    url = "$INSTANCES_URL/$name.json.gz"
    path = download(url)
    curdir = dirname(@__FILE__)
    parentdir = abspath(joinpath(curdir, ".."))
    basedir = abspath(joinpath(parentdir, ".."))
    targetdir = joinpath(basedir, "test/Systems/AC_UC")
    # filename = "$basedir/UC_Cases/$case_name"*"_"*"$date"*".json.gz"
    filename = "$targetdir/$case_name/$date"*".json.gz"
    mkpath(dirname(filename))
    mv(path, filename,force=true)
    return filename
end

function parse_json_case(json_grid_file_path,standard_case_name)

    file = GZip.gzopen(json_grid_file_path)
    json_grid = JSON.parse(file, dicttype = () -> DefaultOrderedDict(nothing))
    GZip.close(file)

    return _create_grid_from_json(json_grid,standard_case_name)
    
end

function _create_grid_from_json(json_grid,case_name)

    grid = PowerGrid(GridID=1)
    grid.S_base = 100
    grid.N_time_steps = json_grid["Parameters"]["Time horizon (h)"]

    for bus in keys(json_grid["Buses"])
        f = replace(bus,"b" => "")
        add_bus!(grid,parse(Int,f),V_max=1.1,V_min=0.9,δ_max=0.6,δ_min=-0.6)
    end

    for load in keys(json_grid["Buses"])
        f = replace(load,"b" => "")
        Pd = first(json_grid["Buses"][load]["Load (MW)"])
        if length(json_grid["Buses"][load]["Load (MW)"]) == 1
            Pd_t = ones(grid.N_time_steps)*json_grid["Buses"][load]["Load (MW)"]
        else
            Pd_t = json_grid["Buses"][load]["Load (MW)"]
        end
        add_load!(grid,parse(Int,f),Pd,0,load_id=parse(Int,f))
        grid.Loads[parse(Int,f)].Pd_t = Pd_t
        grid.Loads[parse(Int,f)].Qd_t = zeros(grid.N_time_steps)
    end

    for gen in keys(json_grid["Generators"])
        gen_id = parse(Int,replace(gen,"g" => ""))
        gen_bus = parse(Int,replace(json_grid["Generators"][gen]["Bus"],"b" => ""))
        C0 = 0
        C1 = sum(json_grid["Generators"][gen]["Production cost curve (\$)"])/length(json_grid["Generators"][gen]["Production cost curve (\$)"])
        C2 = 0
        Pmax = last(json_grid["Generators"][gen]["Production cost curve (MW)"])
        Pmin = json_grid["Generators"][gen]["Production cost curve (MW)"][1]
        C_SU = first(json_grid["Generators"][gen]["Startup costs (\$)"])
        C_SD = 0
        ramp_up = json_grid["Generators"][gen]["Ramp up limit (MW)"]
        ramp_down = json_grid["Generators"][gen]["Ramp down limit (MW)"]
        min_up = json_grid["Generators"][gen]["Minimum uptime (h)"]
        min_down = json_grid["Generators"][gen]["Minimum downtime (h)"]
        Pg = json_grid["Generators"][gen]["Initial power (MW)"]
        Qg = 0
        add_generator!(grid,gen_bus,C0,C1,C2,Pmax,Pmin,0,0,GenType=:base_type,
            start_up_cost=C_SU,shut_down_cost=C_SD,min_up_time=min_up,min_down_time=min_down,Δ_up=ramp_up,
            Δ_down=ramp_down,Pg=Pg,Qg=Qg,gen_id=gen_id)
        grid.Generators[gen_id].Extras = json_grid["Generators"][gen]
    end

    lst = show_cases(false)
    case_id = case_name[2]
    case_name = case_name[1]
    case_name_mod = "pglib_opf_"*case_name
    
    flag = "Normal flow limit (MW)" ∉ keys(json_grid["Transmission lines"]["l1"])

    if flag
        backup_sys = load_system(case_id)
    end
    
    line_count = 0
    for line in keys(json_grid["Transmission lines"])
        line_count += 1
        line_id = parse(Int,replace(line,"l" => ""))
        fbus = parse(Int,replace(json_grid["Transmission lines"][line]["Source bus"],"b"=>""))
        tbus = parse(Int,replace(json_grid["Transmission lines"][line]["Target bus"],"b"=>""))
        r = 0
        x = json_grid["Transmission lines"][line]["Reactance (ohms)"]
        b = 0
        
        if flag
            rating = backup_sys.Branches[line_count].rating
        else
            rating = json_grid["Transmission lines"][line]["Normal flow limit (MW)"]
        end

        add_branch!(grid,fbus,tbus,r_ohms=r,x_ohms=x,b_ohms=b,rating_pu=rating./grid.S_base,branch_id=line_id)
    end

    Y_bus,b_line = Y_Bus_Grid(grid)
    grid.Y_bus = Y_bus
    grid.b_line = b_line

    return grid
end