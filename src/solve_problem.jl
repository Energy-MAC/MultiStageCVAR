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
        SDDP.simulate(model, 1_000, [:pg, :rsv_up, :rsv_dn, :ACE⁺, :ACE⁻])

    path = joinpath(sim_dir, "results")
    open(joinpath(path, "sddp_sol_$(step).json"), "w") do io
        write(io, JSON.json(sims))
    end
    @info "Finished writing SDDP results"
end

function PSI.solve!(problem::PSI.OperationsProblem{MultiStageCVAR})
    if !problem.ext["no_solar"]
        try
            model = problem.ext["mod"]
            SDDP.train(
                model;
                stopping_rules = [SDDP.BoundStalling(20, 1)],
                time_limit = RCVAR.ext["time_limit"],
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
