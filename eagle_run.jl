include("src/MultiStageCVAR.jl")
using CPLEX
using PowerSimulations
using PowerSystems
using Dates

solver = JuMP.optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_MIP_Tolerances_MIPGap" => 2e-6,
    "CPXPARAM_OptimalityTarget" => 1,
)

initial_time = DateTime("2018-04-01T00:00:00")
system_ha = initialize_system("data/HA_sys.json", solver, initial_time)
# system_da = initialize_system("data/DA_sys.json", solver, initial_time)
# system_ed = System("data/RT_sys.json"; time_series_read_only = true)

sddp_solver =
    JuMP.optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Emphasis_Numerical" => 1)

PSY.configure_logging(file_level = Logging.Error, console_level = Logging.Info)

template_hauc = OperationsProblemTemplate(CopperPlatePowerModel)
template_dauc = OperationsProblemTemplate(CopperPlatePowerModel)

for template in [template_dauc, template_hauc]
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    set_device_model!(template, HydroDispatch, FixedOutput)
    set_device_model!(template, ThermalMultiStart, ThermalStandardUnitCommitment)
end

set_service_model!(template_hauc, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
set_service_model!(template_hauc, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))

#DAUC = OperationsProblem(
#    StandardHAUnitCommitmentCC,
#    template_dauc,
#    system_da,
#    optimizer = solver,
#)

#DAUC.ext["cc_restrictions"] = JSON.parsefile("data/cc_restrictions.json")

HAUC = OperationsProblem(
    StandardHAUnitCommitmentCC,
    template_hauc,
    system_ha,
    optimizer = solver,
)

HAUC.ext["cc_restrictions"] = JSON.parsefile("data/cc_restrictions.json")

RCVAR = OperationsProblem(MultiStageCVAR, template_hauc, system_ha, optimizer = sddp_solver)

problems = SimulationProblems(
    #DAUC = DAUC,
    HAUC = HAUC,
    MSCVAR = RCVAR
    #ED = ED,
)

sequence = SimulationSequence(
    problems = problems,
    #feedforward_chronologies = Dict(("DAUC" => "HAUC") => Synchronize(periods = 24)),
    intervals = Dict(
        #"DAUC" => (Hour(24), Consecutive()),
        "HAUC" => (Hour(1), RecedingHorizon()),
        "MSCVAR" => (Hour(1), RecedingHorizon()),
        #    "ED" => (Minute(5), RecedingHorizon()),
    ),
    #feedforward = Dict(
    #    ("ED", :devices, :ThermalMultiStart) => SemiContinuousFF(
    #        binary_source_problem = PSI.ON,
    #        affected_variables = [PSI.ACTIVE_POWER],
    #    ),
    #),
    #ini_cond_chronology = InterProblemChronology(),
    ini_cond_chronology = IntraProblemChronology(),
)

sim = Simulation(
    name = "standard",
    steps = 24,
    problems = problems,
    sequence = sequence,
    initial_time = initial_time,
    simulation_folder = "results",
)

build_out =
    build!(sim; console_level = Logging.Info, file_level = Logging.Error, serialize = false)
execute_out = execute!(sim)
