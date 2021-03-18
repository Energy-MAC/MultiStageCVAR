function PSI.solve!(
    problem::PSI.OperationsProblem{MultiStageCVAR},
)
    model = problem.ext["mod"]
    SDDP.train(model; iteration_limit = 2_000, cut_deletion_minimum = 100)
    sims = SDDP.simulate(
        model,
        1_000,
        [:pg, :rsv_up, :rsv_dn, :ACE⁺, :ACE⁻, :CVAR, :cvar_a]
    )
    # write_model_results!(store, problem, start_time; exports = exports)
    advance_execution_count!(problem)
    return # RunStatus.SUCCESSFUL
end
