# include("DOPF_Constraints.jl")
# A file containting all non-binary/fixed models

function Fixed_UC_Model!(model::Model, grid::PowerGrid, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites, order_book::OrderBook)
    
    unsolved_model = copy(model)
    set_optimizer(unsolved_model, simulation_settings.MILP_solver)

    if length(prerequisites_data.commitable_gen_ids) != 0
        solutions_u_gt = JuMP.value.(model[:u_gt])
        solutions_α_gt = JuMP.value.(model[:α_gt])
        solutions_β_gt = JuMP.value.(model[:β_gt])
        JuMP.fix.(unsolved_model[:u_gt], solutions_u_gt; force = true)
        JuMP.fix.(unsolved_model[:α_gt], solutions_α_gt; force = true)
        JuMP.fix.(unsolved_model[:β_gt], solutions_β_gt; force = true)
        JuMP.unset_binary.(unsolved_model[:u_gt])
        JuMP.unset_binary.(unsolved_model[:α_gt])
        JuMP.unset_binary.(unsolved_model[:β_gt])
    end

    solved_flag_relaxed = optimize_DOPF_model!(unsolved_model, grid, simulation_settings, prerequisites_data, order_book; update_grid=false)
    return unsolved_model, solved_flag_relaxed
end

function Fixed_NCUC_Model!(model::Model, grid::PowerGrid, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites, order_book::OrderBook)
    
    unsolved_model = copy(model)
    set_optimizer(unsolved_model, simulation_settings.MILP_solver)

    if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
        solutions_z_l = JuMP.value.(model[:z_l])
        JuMP.fix.(unsolved_model[:z_l], solutions_z_l; force = true)
        JuMP.unset_binary.(unsolved_model[:z_l])
    end

    if length(prerequisites_data.ac_active_reconf_ids) != 0
        solutions_z_r = JuMP.value.(model[:z_r])
        JuMP.fix.(unsolved_model[:z_r], solutions_z_r; force = true)
        JuMP.unset_binary.(unsolved_model[:z_r])
    end
    
    if length(prerequisites_data.ac_active_coupler_ids) != 0
        solutions_z_c = JuMP.value.(model[:z_c])
        JuMP.fix.(unsolved_model[:z_c], solutions_z_c; force = true)
        JuMP.unset_binary.(unsolved_model[:z_c])
    end

    if length(prerequisites_data.commitable_gen_ids) != 0
        solutions_u_gt = JuMP.value.(model[:u_gt])
        solutions_α_gt = JuMP.value.(model[:α_gt])
        solutions_β_gt = JuMP.value.(model[:β_gt])
        JuMP.fix.(unsolved_model[:u_gt], solutions_u_gt; force = true)
        JuMP.fix.(unsolved_model[:α_gt], solutions_α_gt; force = true)
        JuMP.fix.(unsolved_model[:β_gt], solutions_β_gt; force = true)
        JuMP.unset_binary.(unsolved_model[:u_gt])
        JuMP.unset_binary.(unsolved_model[:α_gt])
        JuMP.unset_binary.(unsolved_model[:β_gt])
    end

    DOPF_single_node_balance!(unsolved_model, grid, simulation_settings, prerequisites_data)
    solved_flag_relaxed = optimize_DOPF_model!(unsolved_model, grid, simulation_settings, prerequisites_data, order_book; update_grid=false)

    return unsolved_model, solved_flag_relaxed
end

function Fixed_SCUC_Model!(model::Model, grid::PowerGrid, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites, order_book::OrderBook)

    unsolved_model = copy(model)
    set_optimizer(unsolved_model, simulation_settings.MILP_solver)

    if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
        solutions_z_l = JuMP.value.(model[:z_l])
        JuMP.fix.(unsolved_model[:z_l], solutions_z_l; force = true)
        JuMP.unset_binary.(unsolved_model[:z_l])
    end

    if length(prerequisites_data.ac_active_reconf_ids) != 0
        solutions_z_r = JuMP.value.(model[:z_r])
        JuMP.fix.(unsolved_model[:z_r], solutions_z_r; force = true)
        JuMP.unset_binary.(unsolved_model[:z_r])
    end
    
    if length(prerequisites_data.ac_active_coupler_ids) != 0
        solutions_z_c = JuMP.value.(model[:z_c])
        JuMP.fix.(unsolved_model[:z_c], solutions_z_c; force = true)
        JuMP.unset_binary.(unsolved_model[:z_c])
    end

    if length(prerequisites_data.commitable_gen_ids) != 0
        solutions_u_gt = JuMP.value.(model[:u_gt])
        solutions_α_gt = JuMP.value.(model[:α_gt])
        solutions_β_gt = JuMP.value.(model[:β_gt])
        JuMP.fix.(unsolved_model[:u_gt], solutions_u_gt; force = true)
        JuMP.fix.(unsolved_model[:α_gt], solutions_α_gt; force = true)
        JuMP.fix.(unsolved_model[:β_gt], solutions_β_gt; force = true)
        JuMP.unset_binary.(unsolved_model[:u_gt])
        JuMP.unset_binary.(unsolved_model[:α_gt])
        JuMP.unset_binary.(unsolved_model[:β_gt])
    end

    if length(keys(prerequisites_data.fixed_commitments)) != 0
        solutions_u_gt_f = JuMP.value.(model[:u_gt_f])
        JuMP.fix.(unsolved_model[:u_gt_f], solutions_u_gt_f; force = true)
    end

    DOPF_single_node_balance!(unsolved_model, grid, simulation_settings, prerequisites_data)
    solved_flag_relaxed = optimize_DOPF_model!(unsolved_model, grid, simulation_settings, prerequisites_data, order_book; update_grid=false)
    return unsolved_model, solved_flag_relaxed
end

function Fixed_SCOPF_Model!(model::Model, grid::PowerGrid, simulation_settings::DOPF_SimulationSettings,
    prerequisites_data::DOPF_Prerequisites, order_book::OrderBook)
    
    unsolved_model = copy(model)
    set_optimizer(unsolved_model, simulation_settings.MILP_solver)
    if length(prerequisites_data.ac_active_dynamic_branch_ids) != 0
        solutions_z_l = JuMP.value.(model[:z_l])
        JuMP.fix.(unsolved_model[:z_l], solutions_z_l; force = true)
        JuMP.unset_binary.(unsolved_model[:z_l])
    end

    if length(prerequisites_data.ac_active_reconf_ids) != 0
        solutions_z_r = JuMP.value.(model[:z_r])
        JuMP.fix.(unsolved_model[:z_r], solutions_z_r; force = true)
        JuMP.unset_binary.(unsolved_model[:z_r])
    end
    
    if length(prerequisites_data.ac_active_coupler_ids) != 0
        solutions_z_c = JuMP.value.(model[:z_c])
        JuMP.fix.(unsolved_model[:z_c], solutions_z_c; force = true)
        JuMP.unset_binary.(unsolved_model[:z_c])
    end

    DOPF_single_node_balance!(unsolved_model, grid, simulation_settings, prerequisites_data)

    solved_flag_relaxed = optimize_DOPF_model!(unsolved_model, grid, simulation_settings, prerequisites_data, order_book; update_grid=false)
    return unsolved_model, solved_flag_relaxed
end

