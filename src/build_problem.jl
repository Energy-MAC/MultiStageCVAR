function apply_cc_constraints!(problem)
    optimization_container = PSI.get_optimization_container(problem)
    restrictions = problem.ext["cc_restrictions"]
    commitment_variables = PSI.get_variable(optimization_container, :On__ThermalMultiStart)
    time_steps = PSI.model_time_steps(optimization_container)
    constraint = PSI.add_cons_container!(
        optimization_container,
        :CC_constraint,
        collect(keys(restrictions)),
        time_steps,
    )
    jump_model = PSI.get_jump_model(optimization_container)
    for t in time_steps, (k, v) in restrictions
        constraint[k, t] =
            JuMP.@constraint(jump_model, sum(commitment_variables[i, t] for i in v) <= 1)
    end
    return
end

function apply_must_run_constraints!(problem)
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    must_run_gens =
        [g for g in PSY.get_components(ThermalMultiStart, system, x -> PSY.get_must_run(x))]
    commitment_variables = PSI.get_variable(optimization_container, :On__ThermalMultiStart)
    for t in time_steps, g in get_name.(must_run_gens)
        JuMP.fix(commitment_variables[g, t], 1.0)
    end
end

function PSI.problem_build!(problem::PSI.OperationsProblem{StandardUnitCommitmentCC})
    PSI.build_impl!(
        PSI.get_optimization_container(problem),
        PSI.get_template(problem),
        PSI.get_system(problem),
    )
    apply_cc_constraints!(problem)
    apply_must_run_constraints!(problem)
end

function set_initial_commitment!(problem::PSI.OperationsProblem{MultiStageCVAR})
    system = PSI.get_system(problem)
    info_map = Dict(
        (PSY.get_name(d) => PSY.get_status(d)) for
        d in PSY.get_components(ThermalMultiStart, system)
    )
    problem.ext["commit_vars"] = info_map
end

function set_initial_reserve!(problem::PSI.OperationsProblem{MultiStageCVAR})
    system = PSI.get_system(problem)
    problem.ext["reserve_vars_up"] = Dict{String, Float64}()
    problem.ext["reserve_vars_dn"] = Dict{String, Float64}()
    for v in problem.ext["balancing_up_devices_names"]
        gen = PSY.get_component(ThermalMultiStart, system, v)
        problem.ext["reserve_vars_up"][v] =
            PSY.get_active_power_limits(gen).max * PSY.get_status(gen)
    end

    for v in problem.ext["balancing_dn_devices_names"]
        gen = PSY.get_component(ThermalMultiStart, system, v)
        problem.ext["reserve_vars_dn"][v] =
            PSY.get_active_power_limits(gen).max * PSY.get_status(gen)
    end
end

function set_time_series!(
    problem::PSI.OperationsProblem{MultiStageCVAR},
    case_initial_time::Dates.DateTime,
)
    @info case_initial_time
    optimization_container = PSI.get_optimization_container(problem)
    time_periods = PSI.model_time_steps(optimization_container)

    problem.ext["total_load"] = zeros(length(time_periods))
    problem.ext["total_solar"] = zeros(length(time_periods))
    problem.ext["total_wind"] = zeros(length(time_periods))
    problem.ext["total_hydro"] = zeros(length(time_periods))

    area = PSY.get_component(Area, system, "1")

    problem.ext["area_solar_forecast_prob"] =
        area_solar_forecast_prob =
            PSY.get_time_series_values(
                Probabilistic,
                area,
                "solar_power";
                start_time = case_initial_time,
            ) / 100

    problem.ext["no_solar"] = all(iszero.(area_solar_forecast_prob))
    problem.ext["no_solar"] && return

    for l in get_components(PowerLoad, system)
        f = get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
        problem.ext["total_load"] .+= f * PSY.get_max_active_power(l)
    end

    for l in
        get_components(RenewableGen, system, x -> get_prime_mover(x) == PrimeMovers.PVe)
        if get_available(l)
            f = get_time_series_values(
                Deterministic,
                l,
                "max_active_power";
                start_time = case_initial_time,
            )
            problem.ext["total_solar"] .+= f * PSY.get_max_active_power(l)
        end
    end

    for l in
        get_components(RenewableGen, system, x -> get_prime_mover(x) != PrimeMovers.PVe)
        if get_available(l)
            f = get_time_series_values(
                Deterministic,
                l,
                "max_active_power";
                start_time = case_initial_time,
            )
            problem.ext["total_wind"] .+= f * PSY.get_max_active_power(l)
        end
    end

    for l in get_components(HydroGen, system)
        f = get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
        problem.ext["total_hydro"] .+= f * PSY.get_max_active_power(l)
    end
end

function set_inputs_dic!(problem::PSI.OperationsProblem{MultiStageCVAR})
    case_initial_time = PSI.get_initial_time(problem)
    problem.ext["MINS_IN_HOUR"] = MINS_IN_HOUR = 60.0
    problem.ext["Δt"] = Δt = 1 / 12
    optimization_container = PSI.get_optimization_container(problem)
    time_periods = PSI.model_time_steps(optimization_container)
    system = get_system(problem)
    thermal_generators = PSY.get_components(PSY.ThermalMultiStart, system)
    problem.ext["thermal_gens_names"] =
        thermal_gens_names = PSY.get_name.(thermal_generators)
    hydro_generators = PSY.get_components(PSY.HydroGen, system)
    hydro_gens_names = PSY.get_name.(hydro_generators)

    problem.ext["cost_component"] = Dict(
        thermal_gens_names .=>
            PSY.get_operation_cost.(thermal_generators) .|> get_variable,
    )

    get_rmp_up_limit(g) = PSY.get_ramp_limits(g).up
    get_rmp_dn_limit(g) = PSY.get_ramp_limits(g).down

    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    reg_up_devices = get_contributing_devices(system, reg_reserve_up)
    reg_dn_devices = get_contributing_devices(system, reg_reserve_dn)
    problem.ext["balancing_up_devices_names"] = get_name.(reg_up_devices)
    problem.ext["balancing_dn_devices_names"] = get_name.(reg_dn_devices)

    problem.ext["pg_lim"] = Dict(
        g => get_active_power_limits(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gens_names
    )

    problem.ext["pg0"] = Dict(
        g => get_active_power(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gens_names
    )

    problem.ext["ramp_up"] = Dict(
        g =>
            get_rmp_up_limit(get_component(ThermalMultiStart, system, g)) *
            Δt *
            MINS_IN_HOUR for g in thermal_gens_names
    )
    problem.ext["ramp_dn"] = Dict(
        g =>
            get_rmp_dn_limit(get_component(ThermalMultiStart, system, g)) *
            Δt *
            MINS_IN_HOUR for g in thermal_gens_names
    )

    problem.ext["not_balancing_thermal"] = intersect(
        thermal_gens_names,
        union(
            problem.ext["balancing_up_devices_names"],
            problem.ext["balancing_dn_devices_names"],
        ),
    )

    set_time_series!(problem, case_initial_time)

    area_solar_forecast_prob = problem.ext["area_solar_forecast_prob"]
    problem.ext["prob_set"] = 1:size(area_solar_forecast_prob)[2]
    problem.ext["prob_values"] =
        (1 / length(problem.ext["prob_set"])) .* ones(size(area_solar_forecast_prob))

    N = length(problem.ext["prob_set"])
    p = fill(1 / N, N)
    problem.ext["MARKOV_TRANSITION"] = Matrix{Float64}[[1.0]', (p ./ sum(p))']
    # Uncomment this code for sparse Transition Matrix
    # ρ = 0.6
    # ev = fill((1 - ρ) / 2, N - 1)
    # ev[1] += ev[1]
    # Tmat = LinearAlgebra.Tridiagonal(reverse(ev), fill(ρ, N), ev)
    Tmat = fill(1 / N, N, N)
    for t in 2:(length(time_periods) - 1)
        push!(problem.ext["MARKOV_TRANSITION"], Tmat)
    end

end

function get_sddp_model(problem::PSI.OperationsProblem{MultiStageCVAR})
    # Checklist
    # TimeSeries generation -> ok
    # Ramp Rates -> ok
    # Device Limits -> ok
    # Probabilistic Forecast -> ok
    # Cost functions -> ok

    settings = PSI.get_settings(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_periods = PSI.model_time_steps(optimization_container)

    thermal_gens_names = [k for (k, v) in problem.ext["commit_vars"] if v]
    balancing_devices_names_up = [k for (k, v) in problem.ext["reserve_vars_up"] if v > 0.0]
    balancing_devices_names_dn = [k for (k, v) in problem.ext["reserve_vars_dn"] if v > 0.0]
    not_balancing_thermal = intersect(
        thermal_gens_names,
        union(balancing_devices_names_up, balancing_devices_names_dn),
    )

    cost_component = problem.ext["cost_component"]
    pg_lim = problem.ext["pg_lim"]
    pg0 = problem.ext["pg0"]

    total_load = problem.ext["total_load"]
    total_solar = problem.ext["total_solar"]
    total_wind = problem.ext["total_wind"]
    total_hydro = problem.ext["total_hydro"]

    area_solar_prob = problem.ext["area_solar_forecast_prob"]

    # There are based on the results of the previous problem run
    commit_vars = problem.ext["commit_vars"]
    reserve_vars_up = problem.ext["reserve_vars_up"]
    reserve_vars_dn = problem.ext["reserve_vars_dn"]

    ramp_up = problem.ext["ramp_up"]
    ramp_dn = problem.ext["ramp_dn"]

    Δt = problem.ext["Δt"]

    model = SDDP.MarkovianPolicyGraph(
        transition_matrices = problem.ext["MARKOV_TRANSITION"],
        sense = :Min,
        lower_bound = 0.0,
        optimizer = sddp_solver,
        direct_mode = true,
    ) do sp, node
        t, markov_state = node[1] - 1, node[2]
        SDDP.set_silent(sp)

        SDDP.@variable(
            sp,
            pg_lim[g].min <= pg[g ∈ thermal_gens_names] <= pg_lim[g].max,
            SDDP.State,
            initial_value = pg0[g]
        )

        SDDP.@variables(
            sp,
            begin
                # Thermal Generation Cost
                cg[thermal_gens_names] >= 0
                # Balacing Reserve variables
                0 <= rsv_up[g ∈ balancing_devices_names_up] <= reserve_vars_up[g]
                0 <= rsv_dn[g ∈ balancing_devices_names_dn] <= reserve_vars_dn[g]
                # ACE
                ACE⁺ >= 0
                ACE⁻ >= 0
                # PWL Cost function auxiliary variables
                0 <= λ[g in thermal_gens_names, i in 1:length(cost_component[g])] <= PSY.get_breakpoint_upperbounds(cost_component[g])[i]
                # slack can be used for debugging purposes. Free if necessary.
                # slack⁺ >= 0
                # slack⁻ >= 0
            end
        )

        SDDP.@constraints(
            sp,
            begin
                # PWL constraints
                [g in thermal_gens_names],
                sum(λ[g, i] for i in 1:length(cost_component[g])) == pg[g].in
                [g in thermal_gens_names],
                cg[g] >= sum(
                    100.0 * λ[g, i] * PSY.get_slopes(cost_component[g])[i] for
                    i in 1:length(cost_component[g])
                )
                # Additional range constraints
                [g in balancing_devices_names_up], pg[g].in + rsv_up[g] <= pg_lim[g].max
                [g in balancing_devices_names_dn], pg[g].in - rsv_dn[g] >= pg_lim[g].min

                [g in balancing_devices_names_up],
                pg[g].out - pg[g].in <= ramp_up[g] - rsv_up[g]
                [g in balancing_devices_names_dn],
                pg[g].in - pg[g].out <= ramp_dn[g] - rsv_dn[g]

                [g in not_balancing_thermal], pg[g].out - pg[g].in <= ramp_up[g]
                [g in not_balancing_thermal], pg[g].in - pg[g].out <= ramp_dn[g]
            end
        )

        ###
        ### Stage-dependent constraints
        ###

        # System balance for the next stage
        if t < length(time_periods)
            SDDP.@constraint(
                sp,
                # total_load[t + 1] + slack⁺ - slack⁻ ==
                0.95*total_load[t + 1] <=
                total_solar[t + 1] +
                total_wind[t + 1] +
                total_hydro[t + 1] +
                #+ slack⁺ - slack⁻ +
                sum(pg[g].out for g in thermal_gens_names)
            )
        else
            # Ignore the generators in the final stage.
            SDDP.@constraint(sp, [g ∈ thermal_gens_names], pg[g].out == pg[g].in)
        end

        ϵ = 0.1
        SDDP.@variable(sp, 0 <= cvar_a, SDDP.State, initial_value = 0.0)
        SDDP.@variable(sp, cvar_y >= 0)
        SDDP.@constraint(sp, cvar_y >= ACE⁺ + ACE⁻ - cvar_a.in)
        SDDP.@expression(sp, CVAR, cvar_a.in + 1 / ϵ * cvar_y)
        if t == 0
            SDDP.@stageobjective(sp, 0.0)
            return
        else
            SDDP.@constraint(
                sp,
                ACE⁺ - ACE⁻ ==
                total_wind[t] +
                total_hydro[t] +
                area_solar_prob[t, markov_state] +
                sum(pg[g].in for g in thermal_gens_names) - total_load[t] +
                sum(rsv_up[g] for g in balancing_devices_names_up) -
                sum(rsv_dn[g] for g in balancing_devices_names_dn)
            )
            SDDP.@stageobjective(sp, sum(cg) * Δt / 1000 +
               1e3 * (ACE⁺ + ACE⁻) +
               sum(rsv_up[g] for g in balancing_devices_names_up) +
               sum(rsv_dn[g] for g in balancing_devices_names_dn)
               )
        end
    end
end

function PSI.problem_build!(problem::PSI.OperationsProblem{MultiStageCVAR})
    set_inputs_dic!(problem)
    set_initial_commitment!(problem)
    set_initial_reserve!(problem)
    mod = get_sddp_model(problem)
    problem.ext["mod"] = mod
end
