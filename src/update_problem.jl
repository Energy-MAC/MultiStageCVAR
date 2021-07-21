function get_jump_var_value(var)
    if JuMP.is_binary(var)
        return round(JuMP.value(var))
    else
        return JuMP.value(var)
    end
end

function update_data!(problem::PSI.OperationsProblem{MultiStageCVAR}, sim::PSI.Simulation)
    current_time = PSI.get_current_time(sim)
    set_time_series!(problem, current_time)
    problem.ext["no_solar"] && return
    @info "Updating SDDP Data"
    problems = PSI.get_problems(sim)
    hauc = problems["HAUC"]
    hauc_optimization_container = PSI.get_optimization_container(hauc)

    variable = PSI.get_variable(hauc_optimization_container, :On__ThermalMultiStart)
    for (k, v) in problem.ext["commit_vars"]
        var = variable[k, 1]
        problem.ext["commit_vars"][k] = get_jump_var_value(var)
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
        problem.ext["reserve_vars_dn"][k] = round(get_jump_var_value(var), digits = 2)
    end

    variable = PSI.get_variable(hauc_optimization_container, :P__ThermalMultiStart)
    for (k, v) in problem.ext["commit_vars"]
        if !v
            problem.ext["pg0"][k] = 0.0
            continue
        end
        var = variable[k, 1]
        var_value = get_jump_var_value(var)
        rev_up = get(problem.ext["reserve_vars_up"], k, 0.0)
        rev_dn = get(problem.ext["reserve_vars_dn"], k, 0.0)
        lims = problem.ext["pg_lim"][k]
        if isapprox(lims.min, var_value - rev_dn)
            problem.ext["pg0"][k] = round(lims.min + rev_dn; digits = 2)
        elseif isapprox(lims.max, var_value + rev_up)
            problem.ext["pg0"][k] = round(lims.max - rev_up; digits = 2)
        else
            problem.ext["pg0"][k] = round(var_value, digits = 2)
        end
        if !(lims.min + rev_dn <= problem.ext["pg0"][k] <= lims.max - rev_up)
            @warn k, lims.min + rev_dn, problem.ext["pg0"][k], lims.max - rev_up
        end
    end
    return
end

function PSI.update_problem!(
    problem::PSI.OperationsProblem{MultiStageCVAR},
    sim::PSI.Simulation,
)
    update_data!(problem, sim)
    if !problem.ext["no_solar"]
        problem.ext["mod"] = get_sddp_model(problem)
    end
end
