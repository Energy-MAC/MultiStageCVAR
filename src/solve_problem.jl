function PSI.write_problem_results!(
    step::Int,
    problem::PSI.OperationsProblem{MultiStageCVAR},
    start_time::Dates.DateTime,
    store::PSI.SimulationStore,
    exports,
)
    problem.ext["no_solar"] && return

    model = problem.ext["mod"]
    sim_dir = dirname(PSI.get_output_dir(problem))
    sims =
        SDDP.simulate(model, 5_000, [:pg, :rsv_up, :rsv_dn, :ACE⁺, :ACE⁻, :CVAR, :cvar_a])

    path = joinpath(sim_dir, "results")
    open(joinpath(path, "sddp_sol_$(step).json"), "w") do io
        write(io, JSON.json(sims))
    end

    # thermal_gens_names = [k for (k, v) in problem.ext["commit_vars"] if v]
    # balancing_devices_names_up = [k for (k, v) in problem.ext["reserve_vars_up"] if v > 0.0]
    # balancing_devices_names_dn = [k for (k, v) in problem.ext["reserve_vars_dn"] if v > 0.0]

    # open(joinpath(path, "sddp_gen_$(step).json"), "w") do io
    #     write(io, JSON.json(thermal_gens_names))
    #     write(io, JSON.json(balancing_devices_names_up))
    #     write(io, JSON.json(balancing_devices_names_dn))
    # end
end

function PSI.solve!(problem::PSI.OperationsProblem{MultiStageCVAR})
    if !problem.ext["no_solar"]
        try
            model = problem.ext["mod"]
            SDDP.train(
                model;
                stopping_rules = [SDDP.BoundStalling(20, 1)],
                time_limit = 10000,
                cut_deletion_minimum = 100,
                run_numerical_stability_report = false,
            )
        catch e
            problem.ext["no_solar"] = true
            @info "step_failed"
        end
    end
    return PSI.RunStatus.SUCCESSFUL
end
