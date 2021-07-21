function read_sddp_results(file::AbstractString)
    res_dict = JSON.parsefile(file)
    return convert(Vector{Vector{Dict{String, Any}}}, res_dict)
end

function initialize_system(system_file::String, solver, initial_time::Dates.DateTime)
    system_ha = System(system_file; time_series_read_only = true)
    interval = get_forecast_interval(system_ha)
    template = OperationsProblemTemplate(CopperPlatePowerModel)
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, HydroDispatch, FixedOutput)
    set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
    set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
    set_device_model!(template, ThermalMultiStart, ThermalBasicUnitCommitment)

    UC =
        OperationsProblem(StandardDAUnitCommitmentCC, template, system_ha, optimizer = solver)
    UC.ext["cc_restrictions"] = JSON.parsefile("data/cc_restrictions.json")

    problems = SimulationProblems(UC = UC)

    sequence = SimulationSequence(
        problems = problems,
        intervals = Dict("UC" => (interval, RecedingHorizon())),
        ini_cond_chronology = IntraProblemChronology(),
    )

    tmp_dir = mktempdir()
    sim = Simulation(
        name = "standard",
        steps = 1,
        problems = problems,
        sequence = sequence,
        initial_time = initial_time,
        simulation_folder = tmp_dir,
    )

    build_out = build!(
        sim;
        console_level = Logging.Info,
        file_level = Logging.Error,
        serialize = false,
    )
    execute_out = execute!(sim)

    results = SimulationResults(joinpath(tmp_dir, "standard"); ignore_status = true)
    res = get_problem_results(results, "UC")
    results = read_realized_variables(res)

    for name in names(results[:On__ThermalMultiStart])
        status = results[:On__ThermalMultiStart][!, name]
        eltype(status) <: DateTime && continue
        th = get_component(ThermalMultiStart, system_ha, name)
        min_lim = get_active_power_limits(th).min
        max_lim = get_active_power_limits(th).max
        set_status!(th, Bool(round(status[1])))
        base_power = get_base_power(th)
        power_ = results[:P__ThermalMultiStart][!, name][1]
        # power_ = status[1]*((min_lim + power[1]))
        power_ = isapprox(power_, max_lim; atol = 1e-3) ? max_lim - 1e-2 : power_
        power_ = isapprox(power_, min_lim; atol = 1e-3) ? min_lim + 1e-3 : power_
        set_active_power!(th, power_)
        set_reactive_power!(th, 0.0)
    end

    to_json(system_ha, system_file; force = true)
    return system_ha
end
