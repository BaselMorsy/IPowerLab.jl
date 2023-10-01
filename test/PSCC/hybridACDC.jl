using IPowerLab
using Ipopt
using Gurobi
using JuMP
using Plots

function get_case(case_id, modularization; load_factor=1, line_factor=1, gen_factor=1, conv_factor=1, load_shedding_cost=1000)
    type=:ACDC_cases
    date="2017-02-18"
    grid = load_system(case_id; type=type, date=date, start_node_from_1=true)
    set_line_capacity_multiplier!(grid, line_factor)
    set_generation_multiplier!(grid, gen_factor)
    for g in keys(grid.Generators)
        if grid.Generators[g].GenType != :virtual
            grid.Generators[g].Pg_min = 0
            grid.Generators[g].Δ_up = 0.2*grid.Generators[g].Pg_max
            grid.Generators[g].Δ_down = grid.Generators[g].Δ_up
        end
    end
    set_load_multiplier!(grid, load_factor)
    set_converter_capacity_multiplier!(grid, conv_factor)
    
    for d in keys(grid.Loads)
        grid.Loads[d].Shedding_Cost = load_shedding_cost
    end
    converter_ids = deepcopy(keys(grid.Converters))
    if modularization != :none
        for conv_id in converter_ids
            split_element!(grid, :Converter, conv_id ;n=2,modularization=modularization)
        end
    else
        modularization = :discrete
    end
    return grid, modularization
end

function solve_case(grid, modularization, DCC, apply_splitting; meta_solver=:none,conv_redispatch_penalty=1000,
    cheat_sheet=[], t_max=100,simulation_type=:SCOPF,relaxed_physics_lines=[], relaxed_physics_nodes=[])
    
    time_horizon = [1]
    meta_solver_instance = MetaSolver()
    meta_solver_instance.misc["t_max"] = t_max
    if cheat_sheet != []
        meta_solver_instance.misc["k_final"] = cheat_sheet["k_final"]
        if "starting_topology" ∈ keys(cheat_sheet)
            meta_solver_instance.misc["starting_topology"] = cheat_sheet["starting_topology"]
        end
    end

    SimulationSettings = DOPF_SimulationSettings(time_horizon = time_horizon,
        ac_grid_model = :Bθ,
        dc_grid_model = :NF, 
        converter_model = :NF_lossless,
        dynamic_converter_control = DCC,
        converter_modularization = modularization,
        converter_redispatch_penalty=conv_redispatch_penalty,
        load_shedding = [:post], # if empty then no load shedding is allowed at all
        transmission_switching = [:post], # if empty then TS will be utilized pre and post contingencies
        substation_switching = Dict("splitting" => [:post], "reconf" => [:pre]),
        activate_temporal_switching = false, # to impose switching frequency constraints
        max_transmission_switching = Dict("pre_contingency" => Inf, "post_contingency" => Inf, "MCDT" => 1),
        max_substation_reconf = Dict("MCDT" => 1),
        max_busbar_splitting = Dict("pre_contingency" => Inf, "post_contingency" => Inf, "MCDT" => 1),

        contingency_types = [:ac_branch,:conv], # all possibilities will be explained in documentation
        contingency_redispatch_condition = :loss_of_generation,

        NLP_solver = Ipopt.Optimizer,
        MILP_solver = Gurobi.Optimizer,
        Meta_solver = meta_solver,
        Meta_solver_instance=meta_solver_instance,
        Parallels = true)

    # Switching specifications
    if apply_splitting
        S = []
        SS = collect(Set([grid.Converters[i].AC_Bus_ID for i in keys(grid.Converters)])) # IDs of buses considered as substations
        for bus_id in SS
            if length(grid.Buses[bus_id].ConnectedLinesIDs) ≥ 2
                push!(S, bus_id)
            end
        end

        for bus_id in keys(grid.Buses)
            if length(grid.Buses[bus_id].ConnectedLinesIDs) ≥ 5 && bus_id ∉ S
                push!(S, bus_id)
            end
        end
    else
        S = []
    end
    H = [] # IDs of buses considered as hybrid substations with b2b Converters
    B2B_cap = Dict() # dictionary (substation_id => capacity) of back to back converter capacity
    ADB = [] # IDs of branches that are allowed to be switched on/off
    FDB = [] # IDs of branches that are switchable but you prefer to give them a specific status based on some knowledge you have
    FT = Dict() # Dictionary (branch_id => status) of fixed topology
    switching_specs = Dict("substation_ids" => S, "hybrid_substation_ids" => H, "B2B_cap" => Dict(),
        "ac_active_dynamic_branch_ids" => ADB, "ac_fixed_dynamic_branch_ids" => FDB, "fixed_topology" => FT)

    # Compile grid topology
    compile_grid_topology!(grid, SimulationSettings, switching_specs)

    all_gen_ids = [gen_id for gen_id in keys(grid.Generators) if grid.Generators[gen_id].GenType != :virtual]
    contingency_specs = Dict("include_leafs" => false, "k" => 1)
    generation_specs = Dict("commitable_gen_ids" => [], "non_commitable_gen_ids" => all_gen_ids, "fixed_commitments" => Dict(), "fixed_schedules" => Dict())
    reference_node = nothing
    prerequisites, virtual_market = compile_prerequisites_DOPF_no_market!(grid, simulation_type, SimulationSettings, contingency_specs, switching_specs,
        generation_specs; relaxed_physics_lines=relaxed_physics_lines, relaxed_physics_nodes=relaxed_physics_nodes, relaxed_capacity_lines=[],reference_node=reference_node)

    solved_flag, model = run_DOPF_simulation_no_market!(grid, SimulationSettings, prerequisites; update_grid=true, return_model=true)
    if meta_solver == :none
        return grid, model, prerequisites
    else
        return grid, SimulationSettings.Meta_solver_instance
    end
    
end



lst = show_cases(true; type=:ACDC_cases)
case_id = 14

line_factor = 1
gen_factor = 2
conv_factor = 0.5

grid, modularization = get_case(case_id, :none; load_factor=1.5, line_factor=line_factor, gen_factor=gen_factor, conv_factor=conv_factor)
solved_grid, meta_solver_cont,prerequisites = solve_case(grid, modularization, true, false; meta_solver=:none,conv_redispatch_penalty=10,
        cheat_sheet=[], t_max=100,simulation_type=:SCOPF, 
        relaxed_physics_lines=collect(keys(grid.Branches)), relaxed_physics_nodes=collect(keys(grid.Buses)))

operating_cost_nomod = []
operating_cost_cont = []
operating_cost_SCOPF = []

solved_grid_nomod_arr = []
solved_grid_cont_arr = []
solved_grid_SCOPF_arr = []

meta_solver_nomod_arr = []
meta_solver_cont_arr = []
models_arr = []

loading_factors = 0.75:0.25:1.5
cheat_sheet_cont = []
cheat_sheet_nomod = []
t_max = 500
meta_solver = :CCG

for loading in loading_factors

    grid_cont, modularization = get_case(case_id, :continuous; load_factor=loading, line_factor=line_factor, gen_factor=gen_factor, conv_factor=conv_factor)
    solved_grid_cont, meta_solver_cont = solve_case(grid_cont, modularization, true, true; meta_solver=meta_solver,conv_redispatch_penalty=10,
        cheat_sheet=cheat_sheet_cont, t_max=t_max, simulation_type=:SCOPF)
    push!(operating_cost_cont, solved_grid_cont.Operating_Cost[1])
    push!(solved_grid_cont_arr, solved_grid_cont)
    push!(meta_solver_cont_arr, meta_solver_cont)
    if meta_solver == :CCG
        cheat_sheet_cont = deepcopy(meta_solver_cont.misc)
    else
        cheat_sheet_cont = []
    end


    grid_nomod, modularization = get_case(case_id, :none; load_factor=loading, line_factor=line_factor, gen_factor=gen_factor, conv_factor=conv_factor)
    if cheat_sheet_nomod == [] && meta_solver == :CCG
        cheat_sheet_nomod = deepcopy(meta_solver_cont.misc)
        delete!(cheat_sheet_nomod,"starting_topology")
    end
    solved_grid_nomod, meta_solver_nomod = solve_case(grid_nomod, modularization, true, true; meta_solver=meta_solver,conv_redispatch_penalty=10,
        cheat_sheet=cheat_sheet_nomod, t_max=t_max,simulation_type=:SCOPF)
    push!(operating_cost_nomod, solved_grid_nomod.Operating_Cost[1])
    push!(solved_grid_nomod_arr, solved_grid_nomod)
    push!(meta_solver_nomod_arr, meta_solver_nomod)
    if meta_solver == :CCG
        cheat_sheet_nomod = deepcopy(meta_solver_nomod.misc)
    else
        cheat_sheet_nomod = []
    end

    grid_SCOPF, modularization = get_case(case_id, :none; load_factor=loading, line_factor=line_factor, gen_factor=gen_factor, conv_factor=conv_factor)
    solved_grid_SCOPF, model, _ = solve_case(grid_SCOPF, modularization, true, false; meta_solver=:none,conv_redispatch_penalty=10,
        cheat_sheet=[], t_max=t_max,simulation_type=:SCOPF)
    push!(operating_cost_SCOPF, solved_grid_SCOPF.Operating_Cost[1])
    push!(solved_grid_SCOPF_arr, solved_grid_SCOPF)
    push!(models_arr, model)
end

order = 1:length(gammas)

p_ls_cont = [maximum([sum(JuMP.value.(meta_solver_cont_arr[i].MP_models[:p_ls_ac])[d,k,1] for d in keys(grid.Loads)) for k in meta_solver_cont_arr[i].misc["k_final"][1]]) for i in order ]
p_ls_nomod = [maximum([sum(JuMP.value.(meta_solver_nomod_arr[i].MP_models[:p_ls_ac])[d,k,1] for d in keys(grid.Loads)) for k in meta_solver_nomod_arr[i].misc["k_final"][1]]) for i in order ]
p_ls_SCOPF = [maximum([sum(JuMP.value.(models_arr[i][:p_ls_ac])[d,k,1] for d in keys(grid.Loads)) for k in prerequisites.k]) for i in order ]


p_curt_cont = [maximum([sum(JuMP.value.(meta_solver_cont_arr[i].MP_models[:p_curt])[g,k,1] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual) for k in meta_solver_cont_arr[i].misc["k_final"][1]]) for i in order ]
p_curt_nomod = [maximum([sum(JuMP.value.(meta_solver_nomod_arr[i].MP_models[:p_curt])[g,k,1] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual) for k in meta_solver_nomod_arr[i].misc["k_final"][1]]) for i in order ]
p_curt_SCOPF = [maximum([sum(JuMP.value.(models_arr[i][:p_curt])[g,k,1] for g in keys(grid.Generators) if grid.Generators[g].GenType != :virtual) for k in prerequisites.k]) for i in order ]


Δp_conv_cont = [maximum([sum(JuMP.value.(meta_solver_cont_arr[i].MP_models[:Δp_conv_pair])[p,k,1] for p in keys(solved_grid_cont_arr[i].Converter_Duplets)) for k in meta_solver_cont_arr[i].misc["k_final"][1]]) for i in order ]
Δp_conv_nomod = [maximum([sum(JuMP.value.(meta_solver_nomod_arr[i].MP_models[:Δp_conv_ac])[g,k,1] for g in keys(solved_grid_nomod_arr[i].Generators) if solved_grid_nomod_arr[i].Generators[g].GenType == :virtual) for k in meta_solver_nomod_arr[i].misc["k_final"][1]]) for i in order ]
Δp_conv_SCOPF = [maximum([sum(JuMP.value.(models_arr[i][:Δp_conv_ac])[g,k,1] for g in keys(grid.Generators) if grid.Generators[g].GenType == :virtual) for k in prerequisites.k]) for i in order ]

Γ_cont = [sum(JuMP.value.(meta_solver_cont_arr[i].MP_models[:Γ])) for i in order]
Γ_nomod = [sum(JuMP.value.(meta_solver_nomod_arr[i].MP_models[:Γ])) for i in order]
Γ_SCOPF = [sum(JuMP.value.(models_arr[i][:Γ])) for i in order]

y = 1
x = 3

plot(gammas[y:x],p_ls_cont[y:x],label="NTR-c", color=:blue)
scatter!(gammas[y:x],p_ls_cont[y:x],color=:blue, label="")
plot!(gammas[y:x],p_ls_nomod[y:x],label="NTR-n", color=:red)
scatter!(gammas[y:x],p_ls_nomod[y:x], color=:red, label="")
plot!(gammas[y:x],p_ls_SCOPF[y:x],label="SCOPF",color=:green)
scatter!(gammas[y:x],p_ls_SCOPF[y:x], color=:green, label="")
xlabel!("Loading level")
ylabel!("Post Contingency Load Shedding (MW)")
savefig("IEEE67_pls.png")


plot(gammas[y:x],p_curt_cont[y:x],label="NTR-c", color=:blue)
scatter!(gammas[y:x],p_curt_cont[y:x],color=:blue, label="")
plot!(gammas[y:x],p_curt_nomod[y:x],label="NTR-n", color=:red)
scatter!(gammas[y:x],p_curt_nomod[y:x],color=:red, label="")
plot!(gammas[y:x],p_curt_SCOPF[y:x],label="SCOPF")
xlabel!("Loading level")
ylabel!("Post Contingency Curtailment (MW)")
savefig("IEEE67_pcurt.png")

plot(gammas[y:x],Δp_conv_cont[y:x],label="NTR-c", color=:blue)
scatter!(gammas[y:x],Δp_conv_cont[y:x],color=:blue, label="")
plot!(gammas[y:x],Δp_conv_nomod[y:x],label="NTR-n", color=:red)
scatter!(gammas[y:x],Δp_conv_nomod[y:x],color=:red, label="")
plot!(gammas[y:x],Δp_conv_SCOPF[y:x],label="SCOPF")
xlabel!("Loading level")
ylabel!("Post Contingency Converter Redispatch (MW)")
savefig("IEEE67_dp_conv.png")


plot(gammas[y:x],operating_cost_cont[y:x]  ,label="NTR-c",color=:blue)
scatter!(gammas[y:x],operating_cost_cont[y:x] ,label="",color=:blue)
plot!(gammas[y:x], operating_cost_nomod[y:x] ,label="NTR-n",color=:red)
scatter!(gammas[y:x], operating_cost_nomod[y:x] ,label="",color=:red)
plot!(gammas[y:x], operating_cost_SCOPF[y:x], label="SCOPF", color=:green)
scatter!(gammas[y:x],operating_cost_SCOPF[y:x],label="",color=:green)
xlabel!("Loading level")
ylabel!("Operating cost")
savefig("IEEE67_opex.png")