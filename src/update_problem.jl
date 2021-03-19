function get_jump_var_value(var)
    if JuMP.is_binary(var)
        return round(JuMP.value(var))
    else
        return JuMP.value(var)
    end
end

function update_data!(problem::PSI.OperationsProblem{MultiStageCVAR}, sim::PSI.Simulation)
    @info "Updating SDDP Data"
    current_time = PSI.get_current_time(sim)
    problems = PSI.get_problems(sim)
    hauc = problems["UC"]
    hauc_optimization_container = PSI.get_optimization_container(hauc)

    variable = PSI.get_variable(hauc_optimization_container, :On__ThermalMultiStart)
    for (k, v) in problem.ext["commit_vars"]
        var = variable[k, 1]

        problem.ext["commit_vars"][k] = get_jump_var_value(var)
    end

    variable = PSI.get_variable(hauc_optimization_container, :P__ThermalMultiStart)
    for (k, v) in problem.ext["pg0"]
        var = variable[k, 1]
        problem.ext["pg0"][k] = get_jump_var_value(var)
    end

    variable =
        PSI.get_variable(hauc_optimization_container, :REG_UP__VariableReserve_ReserveUp)
    for (k, v) in problem.ext["reserve_vars_up"]
        var = variable[k, 1]
        problem.ext["reserve_vars_up"][k] = get_jump_var_value(var)
    end

    variable =
        PSI.get_variable(hauc_optimization_container, :REG_DN__VariableReserve_ReserveDown)
    for (k, v) in problem.ext["reserve_vars_dn"]
        var = variable[k, 1]
        problem.ext["reserve_vars_dn"][k] = get_jump_var_value(var)
    end
    set_time_series!(problem, current_time)
end

function PSI.update_problem!(
    problem::PSI.OperationsProblem{MultiStageCVAR},
    sim::PSI.Simulation,
)
    update_data!(problem, sim)
    problem.ext["mod"] = get_sddp_model(problem::PSI.OperationsProblem{MultiStageCVAR})
end
