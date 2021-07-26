using HDF5
using JSON
using Dates
using PlotlyJS
using PowerSystems
using PowerSimulations
using ColorSchemes
using Colors
using Statistics
const PSI = PowerSimulations
const cscheme =
    vcat(reverse(ColorSchemes.devon.colors[1:5:end]), ColorSchemes.devon.colors[1:5:end])[4:101]
cmap = [[Float64(ix) / 100, "#$(hex(v))"] for (ix, v) in enumerate(cscheme)]
cmap[1][1] = 0.0
cmap[end][1] = 1.0

system_ed = System("data/RT_sys.json"; time_series_read_only = true)

function read_sddp_results(file::AbstractString)
    res_dict = JSON.parsefile(file)
    return convert(Vector{Vector{Dict{String, Any}}}, res_dict)
end

function get_psi_ed_ace(res, system_ed, date, realization_file)
    results = read_realized_variables(res; initial_time = date, len = 24)
    params = read_realized_parameters(res; initial_time = date, len = 24)
    wt_names =
        get_name.(
            get_components(
                RenewableGen,
                system_ed,
                x -> get_prime_mover(x) == PrimeMovers.WT,
            ),
        )

    total_thermal = zeros(24)
    for c in eachcol(results[:P__ThermalMultiStart])
        eltype(c) != Float64 && continue
        total_thermal .+= c
    end

    total_load = zeros(24)
    for c in eachcol(params[:P__max_active_power__PowerLoad])
        eltype(c) != Float64 && continue
        total_load .+= c
    end

    total_hydro = zeros(24)
    for c in eachcol(params[:P__max_active_power__HydroDispatch])
        eltype(c) != Float64 && continue
        total_hydro .+= c
    end

    total_wind = zeros(24)
    for g in wt_names
        total_wind .+= params[:P__max_active_power__RenewableDispatch][!, g]
    end

    day_count = dayofyear(date) - 1
    hour = Dates.hour(date)
    ix = (day_count * 24 + hour) * 12
    actual = read(realization_file, "Power")
    @show actual_sample = actual[(ix + 1):(ix + 24)] / 100

    net_power =
        -1 .* total_load .- total_hydro .- total_wind .- total_thermal .- actual_sample

    total_res_up = zeros(24)
    for c in eachcol(results[:REG_UP__VariableReserve_ReserveUp])
        eltype(c) != Float64 && continue
        total_res_up .+= c
    end

    total_res_down = zeros(24)
    for c in eachcol(results[:REG_DN__VariableReserve_ReserveDown])
        eltype(c) != Float64 && continue
        total_res_down .+= c
    end

    ACE = zeros(24)
    for (ix, val) in enumerate(net_power)
        if val < 0.0
            ACE[ix] = min(0.0, val + total_res_up[ix]) * 0.1
        elseif val > 0.00
            ACE[ix] = -1 * min(0.0, total_res_down[ix] - val) * 0.1
        end
    end
    return ACE
end

function make_ace_plots(file, date, res, system_ed, realization_file, color_bar = false)
    sddp_results = read_sddp_results(file)

    mat = Matrix{Float64}(undef, 24, 1000)
    for t in 1:24, s in 1:1000
        data = sddp_results[s][t]
        mat[t, s] = data["ACE⁺"] - data["ACE⁻"]
    end

    qmat = Matrix{Float64}(undef, 24, 99)
    for t in 1:24
        t_ = t + 1 > 24 ? 24 : t + 1
        qmat[t, :] = Statistics.quantile(mat[t_, :], 0.01:0.01:0.99)
    end
    date_range = range(date; length = 24, step = Minute(5))
    traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
    for i in 1:1:98
        trace = scatter(;
            x = vcat(date_range, reverse(date_range)),
            y = vcat(qmat[:, i], reverse(qmat[:, i + 1])) ./ 10,
            fill = "tonextx",
            showlegend = false,
            mode = "none",
            fillcolor = cscheme[i],
            smoothing = 1.3,
            #opacity  = (50-i)/50,
            #marker=attr(color="black", line_width=0, size=4, symbol="circle")
        )
        push!(traces, trace)
    end
    if color_bar
        colorbar_trace = scatter(
            x = [0.0],
            y = [0.0],
            showlegend = false,
            marker = attr(
                colorbar = attr(
                    tickmode = "array",
                    tickvals = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6],
                    ticktext = ["20th", "30th", "40th", "50th", "60th", "70th", "80th"],
                    title = attr(text = "Percentile", font = attr(size = 18)),
                    titleside = "right",
                    tickfont = attr(size = 18),
                ),
                colorscale = cmap,
            ),
        )
        push!(traces, colorbar_trace)
    end

    ub_trace = scatter(;
        x = date_range,
        y = qmat[:, 20] .* 0.1,
        name = "Bottom 20th Percentile",
        model = Line,
        line_color = "#ff7f0e",
    )
    push!(traces, ub_trace)
    lb_trace = scatter(;
        x = date_range,
        y = qmat[:, 80] .* 0.1,
        name = "Top 20th Percentile",
        model = Line,
        line_color = "#d62728",
    )
    push!(traces, lb_trace)
    layout = Layout(;
        # title = "Total Area Probabilistic Forecast",
        xaxis_range = [date_range[1], date_range[end]],
        xaxis = attr(
            tickfont = attr(size = 16),
            title = attr(text = "Forecast Timestamp", standoff = 0, font = attr(size = 18)),
        ),
        yaxis = attr(
            range = [-3, 3],
            title = attr(text = "ACE [GW]", font = attr(size = 18)),
            tickfont = attr(size = 16),
            #autorange = true,
            showgrid = true,
            nticks = 7,
            zeroline = true,
            font = attr(size = 18),
        ),
        legend = attr(orientation = "h", x = 0, xanchor = "left", y = 1),
        # coloraxis= attr(colorbar = attr(title = "Percentile"), colorscale = cmap)
    )
    p = plot(traces, layout)
    return p
end

results = SimulationResults("results/standard/", 2; ignore_status = true)
res = get_problem_results(results, "ED")

morning_file = "/Users/jdlara/Dropbox/Code/MultiStageCVAR/results/1/sddp_sol_8.json"
midday_file = "/Users/jdlara/Dropbox/Code/MultiStageCVAR/results/1/sddp_sol_13.json"
afternoon_file = "/Users/jdlara/Dropbox/Code/MultiStageCVAR/results/1/sddp_sol_20.json"
realization_file = h5open(
    "/Users/jdlara/cache/blue_texas/input_data/Solar/5-minute solar actuals/ERCOT132.h5",
)

date_morning = DateTime("2018-04-01T07:00:00")
date_midday = DateTime("2018-04-01T12:00:00")
date_afternoon = DateTime("2018-04-01T19:00:00")

# system_ed = System("data/RT_sys.json"; time_series_read_only = true)
p1 = make_ace_plots(morning_file, date_morning, res, system_ed, realization_file, false)
p2_ace = make_ace_plots(midday_file, date_midday, res, system_ed, realization_file, true)
p3 = make_ace_plots(afternoon_file, date_afternoon, res, system_ed, realization_file, false)

[p1, p2_ace, p3]

savefig([p1, p2_ace, p3], "ACE_output.pdf"; width = 900, height = 1800)

function make_reserve_plots(file, date, res, system_ed, realization_file, color_bar = true)
    sddp_results = read_sddp_results(file)

    mat = Matrix(undef, 24, 1000)
    for t in 1:24, s in 1:1000
        data = sddp_results[s][t]
        mat[t, s] = sum(data["rsv_up"]) - sum(data["rsv_dn"])
    end

    #=
    results = read_realized_variables(res; initial_time = date, len = 24)
    params = read_realized_parameters(res; initial_time = date, len = 24)
    wt_names = get_name.(get_components(RenewableGen, system_ed, x-> get_prime_mover(x) == PrimeMovers.WT))

    total_thermal = zeros(24)
    for c in eachcol(results[:P__ThermalMultiStart])
        eltype(c) != Float64 && continue
        total_thermal .+= c
    end

    total_load = zeros(24)
    for c in eachcol(params[:P__max_active_power__PowerLoad])
        eltype(c) != Float64 && continue
        total_load .+= c
    end

    total_hydro = zeros(24)
    for c in eachcol(params[:P__max_active_power__HydroDispatch])
        eltype(c) != Float64 && continue
        total_hydro .+= c
    end

    total_wind = zeros(24)
    for g in wt_names
        total_wind .+= params[:P__max_active_power__RenewableDispatch][!, g]
    end

    day_count = dayofyear(date) - 1
    hour = Dates.hour(date)
    ix = (day_count*24 + hour )*12
    @show DateTime("2018-01-01") + ix * Minute(5)
    actual = read(realization_file, "Power")
    actual_sample = actual[ix + 1:ix + 24]/100

    @show net_power = -1 .* total_load .- total_hydro .- total_wind .- total_thermal .- actual_sample

    total_res_up = zeros(24)
    for c in eachcol(results[:REG_UP__VariableReserve_ReserveUp])
        eltype(c) != Float64 && continue
        total_res_up .+= c
    end

    total_res_down = zeros(24)
    for c in eachcol(results[:REG_DN__VariableReserve_ReserveDown])
        eltype(c) != Float64 && continue
        total_res_down .+= c
    end

    res_dep = zeros(24)
    for (ix, val) in enumerate(net_power)
        if isapprox(0.0, val, atol= 1e-3)
            print(ix)
            continue
        end
        if val < 0.0
            res_dep[ix] = max(total_res_up[ix], val - total_res_up[ix])
        elseif val > 0.00
            res_dep[ix] = -1*min(total_res_down[ix], val + total_res_down[ix])
        end
    end
    =#
    qmat = Matrix{Float64}(undef, 24, 99)
    for t in 1:24
        t_ = t + 1 > 24 ? 24 : t + 1
        qmat[t, :] = Statistics.quantile(mat[t_, :], 0.01:0.01:0.99)
    end
    date_range = range(date; length = 24, step = Minute(5))
    traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
    for i in 1:1:98
        trace = scatter(;
            x = vcat(date_range, reverse(date_range)),
            y = vcat(qmat[:, i], reverse(qmat[:, i + 1])) ./ 10,
            fill = "tonextx",
            showlegend = false,
            mode = "none",
            fillcolor = cscheme[i],
            smoothing = 1.3,
            #opacity  = (50-i)/50,
            #marker=attr(color="black", line_width=0, size=4, symbol="circle")
        )
        push!(traces, trace)
    end
    if color_bar
        colorbar_trace = scatter(
            x = date,
            y = [0.0],
            showlegend = false,
            marker = attr(
                colorbar = attr(
                    tickmode = "array",
                    tickvals = [1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6],
                    ticktext = ["20th", "30th", "40th", "50th", "60th", "70th", "80th"],
                    title = attr(text = "Percentile", font = attr(size = 18)),
                    titleside = "right",
                    tickfont = attr(size = 18),
                ),
                colorscale = cmap,
            ),
        )
        push!(traces, colorbar_trace)
    end
    ub_trace = scatter(;
        x = date_range,
        y = qmat[:, 25] .* 0.1,
        name = "Bottom 20th Percentile",
        model = Line,
        line_color = "#ff7f0e",
    )
    push!(traces, ub_trace)
    lb_trace = scatter(;
        x = date_range,
        y = qmat[:, 75] .* 0.1,
        name = "Top 20th Percentile",
        model = Line,
        line_color = "#d62728",
    )
    push!(traces, lb_trace)
    layout = Layout(;
        # title = "Total Area Probabilistic Forecast",
        xaxis_range = [date_range[1], date_range[end]],
        xaxis = attr(
            tickfont = attr(size = 16),
            title = attr(text = "Forecast Timestamp", standoff = 0, font = attr(size = 18)),
        ),
        yaxis = attr(
            title = attr(text = " Reserve Deployment [GW]", font = attr(size = 18)),
            tickfont = attr(size = 16),
            range = [-0.7, 0.7],
            # autorange = true,
            nticks = 7,
            showgrid = true,
            zeroline = true,
            font = attr(size = 18),
        ),
        legend = attr(orientation = "h", x = 0, xanchor = "left", y = 1),
        # coloraxis= attr(colorbar = attr(title = "Percentile"), colorscale = cmap)
    )
    p = plot(traces, layout)
end

date_morning = DateTime("2018-04-01T07:00:00")
date_midday = DateTime("2018-04-01T12:00:00")
date_afternoon = DateTime("2018-04-01T19:00:00")

results = SimulationResults("results/standard/", 2; ignore_status = true)
res = get_problem_results(results, "ED")
# system_ed = System("data/RT_sys.json"; time_series_read_only = true)

p1 = make_reserve_plots(morning_file, date_morning, res, system_ed, realization_file, false)

p2 = make_reserve_plots(midday_file, date_midday, res, system_ed, realization_file, true)
display(p2)

p3 = make_reserve_plots(
    afternoon_file,
    date_afternoon,
    res,
    system_ed,
    realization_file,
    false,
)

[p1, p2, p3]

savefig([p1, p2, p3], "rsv_output.pdf"; width = 900, height = 1800)

savefig([p2_ace, p2], "corrected_rsv.pdf"; width = 900, height = 1012)

[p2_ace, p2]

sddp_results = read_sddp_results(afternoon_file);

mat = Matrix(undef, 24, 1000)
for t in 1:24, s in 1:1000
    data = sddp_results[s][t]
    mat[t, s] = sum(v["out"] for v in values(data["pg"]))
end
date = date_afternoon

results = read_realized_variables(res; initial_time = date, len = 24)
params = read_realized_parameters(res; initial_time = date, len = 24)
wt_names =
    get_name.(
        get_components(RenewableGen, system_ed, x -> get_prime_mover(x) == PrimeMovers.WT),
    )

total_thermal = zeros(24)
for c in eachcol(results[:P__ThermalMultiStart])
    eltype(c) != Float64 && continue
    total_thermal .+= c
end

traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
push!(
    traces,
    scatter(;
        x = [date],
        y = [0.00],
        name = "Probabilistic Reserve Deployment",
        mode = "markers",
        marker = attr(
            color = "#1f77b4",
            opacity = 0.5,
            size = 6,
            symbol = "square",
            line = attr(color = "#1f77b4", width = 0.0),
        ),
    ),
)
for i in 1:24
    trace = box(; #x = 1:24,
        y = mat[i, :] * 100,
        mode = "markers",
        marker_color = "#1f77b4",
        marker_size = 2,
        line_width = 1,
        name = date + (i - 1) * Minute(5),
        # boxpoints ="Outliers",
        showlegend = false,
        fillcolor = "cls",
        #marker=attr(color="black", line_width=0, size=4, symbol="circle")
    )
    push!(traces, trace)
end
realization_trace = scatter(;
    x = range(date; length = 24, step = Minute(5)),
    y = total_thermal * 100,
    name = "realization",
    model = Line,
    line_color = "#ff7f0e",
)

push!(traces, realization_trace)
layout = Layout(;
    # title = "Total Area Probabilistic Forecast",
    xaxis = attr(
        title = attr(text = "Forecast Timestamp", standoff = 0, font = attr(size = 18)),
        tickfont = attr(size = 16),
    ),
    yaxis = attr(
        title = "Reserve Deployment [MW]",
        tickfont = attr(size = 16),
        autorange = true,
        showgrid = true,
        zeroline = true,
        font = attr(size = 18),
    ),
    legend = attr(orientation = "h", x = 0, xanchor = "left", y = 1),
)
p = plot(traces, layout)

rec = list_simulation_events(
    PSI.ProblemExecutionEvent,
    "./results/standard/1",
    x -> x.problem == 2,
    wall_time = true,
)
comp_stamp = [r.common.timestamp for r in rec]
sddp_sol_time = Vector()
for i in 1:2:47
    delta_t = comp_stamp[i + 1] - comp_stamp[i]
    push!(sddp_sol_time, delta_t.value)
end

t_range = range(DateTime("2018-01-01"); step = Hour(1), length = 24)[sddp_sol_time .> 1000]

layout_time = Layout(;
    # title = "Total Area Probabilistic Forecast",
    xaxis = attr(
        title = attr(text = "Timestamp", standoff = 0, font = attr(size = 18)),
        tickfont = attr(size = 16),
    ),
    yaxis = attr(
        title = attr(text = "Solution Time [x1000 s]", font = attr(size = 16)),
        autorange = true,
        showgrid = true,
        zeroline = true,
        tickfont = attr(size = 16),
    ),
    legend = attr(orientation = "h", x = 0, xanchor = "left", y = 1),
)

p = plot(
    bar(; x = t_range, y = sddp_sol_time[sddp_sol_time .> 1000] / 1000000),
    layout_time,
)

savefig(p, "sddp_times.pdf"; width = 900, height = 400)

plot(bar(; x = collect()y = sddp_sol_time[sddp_sol_time .> 1000] / 1000))

trace1 = scatter(;
    x = vcat(1:10, 10:-1:1),
    y = vcat(2:11, 9:-1:0),
    fill = "tozerox",
    fillcolor = "rgba(0, 100, 80, 0.2)",
    line_color = "transparent",
    name = "Fair",
    showlegend = false,
)

plot(trace1)




using PowerSystems
using PlotlyJS
using Dates
using HDF5
using Statistics
using ColorSchemes
using Colors

prob_file = h5open("/Users/jdlara/cache/blue_texas/input_data/Solar/ERCOT132.h5")
actual_file = h5open(
    "/Users/jdlara/cache/blue_texas/input_data/Solar/5-minute solar actuals/ERCOT132.h5",
)

dates = [
    DateTime("2018-04-01T07:00:00"),
    DateTime("2018-04-01T12:00:00"),
    DateTime("2018-04-01T19:00:00"),
]

plots = Vector()
cscheme =
    vcat(reverse(ColorSchemes.devon.colors[1:5:end]), ColorSchemes.devon.colors[1:5:end])[4:101]
cmap = [[Float64(ix) / 100, "#$(hex(v))"] for (ix, v) in enumerate(cscheme)]
cmap[1][1] = 0.0
cmap[end][1] = 1.0
colorbar = true
for date in dates
    day_count = dayofyear(date) - 1
    hour = Dates.hour(date)
    ix = (day_count * 24 + hour) * 12
    DateTime("2018-01-01") + ix * Minute(5)
    data_prob = read(prob_file, "Power")
    forecast = data_prob[ix + 1, :, :]
    actual = read(actual_file, "Power")
    actual_sample = actual[(ix + 1):(ix + 24)]
    date_range = range(date; length = 24, step = Minute(5))
    traces = Vector{GenericTrace{Dict{Symbol, Any}}}()
    for i in 1:1:98
        trace = scatter(;
            x = vcat(date_range, reverse(date_range)),
            y = vcat(forecast[:, i] .- 0.01, reverse(forecast[:, i + 1]) .+ 0.01) ./ 1000,
            fill = "tonextx",
            showlegend = false,
            mode = "none",
            fillcolor = cscheme[i],
            smoothing = 1.3,
            yaxis_range = [date_range[1], date_range[end]],
        )
        push!(traces, trace)
    end
    if colorbar
        lb_trace = scatter(;
            x = date_range,
            y = forecast[:, 80] ./ 1000,
            model = Line,
            line_width = 0.00,
            showlegend = false,
            line_color = "black",
        )
        push!(traces, lb_trace)
        colorbar_trace = scatter(
            x = [0.0],
            y = [0.0],
            showlegend = false,
            marker = attr(
                colorbar = attr(
                    tickmode = "array",
                    tickvals = [2.2, 2.74, 3.4, 4, 4.6, 5.14, 5.8],
                    ticktext = ["20th", "30th", "40th", "50th", "60th", "70th", "80th"],
                    title = attr(text = "Percentile", font = attr(size = 18)),
                    titleside = "right",
                    tickfont = attr(size = 18),
                ),
                colorscale = cmap,
            ),
        )
        push!(traces, colorbar_trace)
        colorbar = false
    end

    l = Layout(;
        # title = "Total Area Probabilistic Forecast",
        xaxis = attr(
            tickfont = attr(size = 16),
            title = attr(text = "Forecast Timestamp", standoff = 0, font = attr(size = 18)),
        ),
        yaxis = attr(
            title = attr(text = "Solar Power [GW]", font = attr(size = 18)),
            tickfont = attr(size = 16),
            autorange = true,
            showgrid = true,
            zeroline = true,
            font = attr(size = 18),
        ),
        ygap = 0,
        width = 300,
        height = 550,
        # legend = attr(orientation = "h", x = 0, xanchor = "left", y = 1),
        xaxis_range = [date_range[1], date_range[end]],
    )
    push!(plots, plot(traces, l))
end
[plots...]

savefig([plots...], "output_filename.pdf"; width = 900, height = 1600)
