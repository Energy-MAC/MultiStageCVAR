include("src/MultiStageCVAR.jl")
using CPLEX
using PowerSimulations
using PowerSystems
using Dates
using PlotlyJS

solver = JuMP.optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_MIP_Tolerances_MIPGap" => 1e-2,
    "CPXPARAM_Emphasis_MIP" => 1,
)

initial_time = DateTime("2018-04-01T00:00:00")
system_ha = initialize_system("data/HA_sys.json", solver, initial_time)
system_da = initialize_system("data/DA_sys.json", solver, initial_time)
#system_da = initialize_system("data/DA_sys.json", solver, initial_time)
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
    set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
    set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
end

set_device_model!(template_hauc, ThermalMultiStart, ThermalBasicUnitCommitment)

DAUC = OperationsProblem(
    StandardDAUnitCommitmentCC,
    template_dauc,
    system_da,
    optimizer = solver,
    optimizer_log_print = true,
    warm_start = false
)

DAUC.ext["cc_restrictions"] = JSON.parsefile("data/cc_restrictions.json")

#=
build!(DAUC; output_dir = mktempdir(), serialize = false)
solve!(DAUC)
results = ProblemResults(DAUC)

reg = results.variable_values[:REG_UP__VariableReserve_ReserveUp]
traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
for i in eachcol(reg2)
    if eltype(i) == Float64 && sum(i) > 1e-3
        push!(traces, PlotlyJS.scatter(y = i))
    end
end

plot(traces, Layout())
=#

HAUC = OperationsProblem(
    StandardHAUnitCommitmentCC,
    template_hauc,
    system_ha,
    optimizer = solver,
    optimizer_log_print = false,
    services_slack_variables = true,
    warm_start = false
)

HAUC.ext["cc_restrictions"] = JSON.parsefile("data/cc_restrictions.json")

RCVAR = OperationsProblem(MultiStageCVAR, template_hauc, system_ha, optimizer = sddp_solver)

problems = SimulationProblems(
    #DAUC = DAUC,
    HAUC = HAUC,
    #MSCVAR = RCVAR
    #ED = ED,
)

sequence = SimulationSequence(
    problems = problems,
    #feedforward_chronologies = Dict(("DAUC" => "HAUC") => Synchronize(periods = 24)),
    intervals = Dict(
        "DAUC" => (Hour(24), Consecutive()),
        "HAUC" => (Hour(1), RecedingHorizon()),
        #"MSCVAR" => (Hour(1), RecedingHorizon()),
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
    steps = 1,
    problems = problems,
    sequence = sequence,
    initial_time = initial_time,
    simulation_folder = "results",
)

build_out =
    build!(sim; console_level = Logging.Info, file_level = Logging.Error, serialize = false)
execute_out = execute!(sim)

#=

build!(DAUC; output_dir = mktempdir(), serialize = false)
solve!(DAUC)
results = ProblemResults(DAUC)

reg = results.variable_values[:REG_UP__VariableReserve_ReserveUp]
traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
for i in eachcol(reg)
    if eltype(i) == Float64 && sum(i) > 1e-3
        push!(traces, scatter(y = i))
    end
end
plot(traces, Layout())
=#
