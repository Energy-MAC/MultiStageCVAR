function PSI.write_problem_results!(
    step::Int,
    problem::PSI.OperationsProblem{MultiStageCVAR},
    start_time::Dates.DateTime,
    store::PSI.SimulationStore,
    exports,
)
    model = problem.ext["mod"]
    sim_dir = dirname(PSI.get_output_dir(problem))
    sims =
        SDDP.simulate(model, 1_000, [:pg, :rsv_up, :rsv_dn, :ACE⁺, :ACE⁻, :CVAR, :cvar_a])

    path = joinpath(sim_dir, "results")
    open(joinpath(path, "sddp_sol_$(step).json"), "w") do io
        write(io, JSON.json(sims))
    end
end

function PSI.solve!(problem::PSI.OperationsProblem{MultiStageCVAR})
    model = problem.ext["mod"]
    SDDP.train(
        model;
        print_level = 2,
        stopping_rules = [SDDP.BoundStalling(5, 10)],
        time_limit = 1000,
        cut_deletion_minimum = 100,
    )

    return PSI.RunStatus.SUCCESSFUL
end
