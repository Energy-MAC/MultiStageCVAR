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

    variable =
        PSI.get_variable(hauc_optimization_container, :SPIN__VariableReserve_ReserveUp)
    for (k, v) in problem.ext["reserve_spin"]
        var = variable[k, 1]
        problem.ext["reserve_spin"][k] = round(get_jump_var_value(var), digits = 2)
    end

    variable = PSI.get_variable(hauc_optimization_container, :P__ThermalMultiStart)
    total_p0 = 0.0
    total_p0_ = 0.0
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
        total_p0_ += var_value
        problem.ext["pg0"][k] = clamp(var_value, lims.min, lims.max)
        if !(lims.min + rev_dn <= problem.ext["pg0"][k] <= lims.max - rev_up)
            @warn k, lims.min + rev_dn, problem.ext["pg0"][k], lims.max - rev_up
        end
        total_p0 += problem.ext["pg0"][k]
    end
    @info "init diff" total_p0 - problem.ext["nl_t0"]
    @info "init diff _" total_p0_ - problem.ext["nl_t0"]

    hour = string(problem.ext["hour"] - 1)
    @show hour_string = length(hour) > 1 ? hour : "0" * hour
    mat = CSV.File("/Users/jdlara/Dropbox/Code/MultiStageCVAR/data/tm_online_08162018_$(hour_string)00.csv") |> DataFrame |> Matrix
    Tmat = mat[:, 2:end] .+ 1e-3
    for i in 1:99
        Tmat[i, :] ./= sum(Tmat[i, :])
    end
    steps = PSI.model_time_steps(PSI.get_optimization_container(problems["MSCVAR"]))
    for ix in 2:steps[end]
        problem.ext["MARKOV_TRANSITION"][ix] = Tmat
    end
    return
end

function PSI.update_problem!(
    problem::PSI.OperationsProblem{MultiStageCVAR},
    sim::PSI.Simulation,
)
    problem.ext["hour"] = Hour(PSI.get_current_time(sim)).value
    update_data!(problem, sim)
    if !problem.ext["no_solar"]
        problem.ext["mod"] = get_sddp_model(problem)
    end
end

function PSI.update_problem!(
    problem::PSI.OperationsProblem{StandardHAUnitCommitmentCC},
    sim::PSI.Simulation,
)
    PSI._update_problem!(problem, sim)
    optimization_container = PSI.get_optimization_container(problem)
    res_dn_var =
        PSI.get_variable(optimization_container, :REG_DN__VariableReserve_ReserveDown)
    res_up_var =
        PSI.get_variable(optimization_container, :REG_UP__VariableReserve_ReserveUp)
    spi = PSI.get_variable(optimization_container, :SPIN__VariableReserve_ReserveUp)
    time_steps = PSI.model_time_steps(optimization_container)
    map_v = Dict(:reg_up_da => res_up_var, :reg_dn_da => res_dn_var, :spin_da => spi)
    current_time = PSI.get_current_time(sim)
    simulation_hour = Hour(current_time).value
    for (key, var) in map_v
        for name in axes(var)[1]
            data_ = problem.ext["resv_dauc"][key][!, Symbol(name)]
            for t in time_steps
                if t < 13
                    data = data_[simulation_hour + 1]
                elseif t >= 12
                    data = data_[simulation_hour + 2]
                end
                if data > 1e-3
                    JuMP.set_upper_bound(var[name, t], data * 2)
                    JuMP.set_lower_bound(var[name, t], data)
                elseif isapprox(data, 0.0, atol = 1e-3)
                    JuMP.set_lower_bound(var[name, t], 0.0)
                    JuMP.set_upper_bound(var[name, t], 10.0)
                else
                    error("bad reserve read")
                end
                @assert JuMP.lower_bound(var[name, t]) <= JuMP.upper_bound(var[name, t])
            end
        end
    end
end
