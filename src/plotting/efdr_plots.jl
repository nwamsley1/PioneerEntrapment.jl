using Plots
using DataFrames
using Statistics

# Try to configure GR backend for vector-friendly PDF text (editable in Illustrator)
try
    import GR
    # 0 = string (stroke), 1 = character outline, 2 = filled, 3 = polygon
    GR.setcharquality(0)
catch
    # If GR is not the active backend, this is a no-op
end

function _auto_axis_limits(df::DataFrame, qval_col::Symbol, efdr_cols::Vector{Symbol})
    xvals = skipmissing(df[!, qval_col])
    xmax = maximum(xvals; init=0.0)
    ymax = 0.0
    for c in efdr_cols
        ymax = max(ymax, maximum(skipmissing(df[!, c]); init=0.0))
    end
    # pad slightly, and clip to a reasonable upper bound
    xupper = min(max(0.01, 1.05 * xmax), 0.1)
    yupper = min(max(0.01, 1.05 * ymax), 0.1)
    return (0.0, xupper), (0.0, yupper)
end

function plot_efdr_vs_qval(df::DataFrame, qval_col::Symbol, efdr_cols::Vector{Symbol};
                          title="Entrapment vs Decoy FDR",
                          xlabel="Decoy FDR",
                          ylabel="Entrapment FDR",
                          labels=nothing,
                          colors=nothing,
                          xlims=nothing,
                          ylims=nothing,
                          legend=:bottomright,
                          diagonal=true,
                          linewidth::Real=1.5)
    sorted_indices = sortperm(df[!, qval_col])
    sorted_df = df[sorted_indices, :]
    if isnothing(xlims) || isnothing(ylims)
        ax_x, ax_y = _auto_axis_limits(sorted_df, qval_col, efdr_cols)
        xlims = isnothing(xlims) ? ax_x : xlims
        ylims = isnothing(ylims) ? ax_y : ylims
    end
    p = plot(title=title, xlabel=xlabel, ylabel=ylabel, xlims=xlims, ylims=ylims, legend=legend, size=(600, 500), dpi=300, fontfamily="Helvetica")
    if diagonal
        max_val = min(xlims[2], ylims[2])
        plot!(p, [0, max_val], [0, max_val], label="y=x", linestyle=:dash, color=:gray, alpha=0.5)
    end
    if isnothing(colors)
        colors = [:blue, :red, :green, :orange, :purple]
    end
    for (i, efdr_col) in enumerate(efdr_cols)
        label = isnothing(labels) ? String(efdr_col) : labels[i]
        color = colors[mod1(i, length(colors))]
        plot!(p, sorted_df[!, qval_col], sorted_df[!, efdr_col], label=label, color=color, linewidth=linewidth, alpha=0.9)
    end
    return p
end

function plot_efdr_comparison(df::DataFrame, score_col::Symbol, qval_col::Symbol; method_types::Vector=[CombinedEFDR, PairedEFDR], kwargs...)
    efdr_cols = Symbol[]
    labels = String[]
    for method_type in method_types
        method_name = method_type == CombinedEFDR ? "combined" : method_type == PairedEFDR ? "paired" : error("Unknown EFDR method type: $method_type")
        efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
        if !hasproperty(df, efdr_col)
            error("Column $efdr_col not found. Make sure to run add_efdr_columns! first.")
        end
        push!(efdr_cols, efdr_col)
        push!(labels, titlecase(method_name) * " EFDR")
    end
    default_kwargs = Dict(:title => "EFDR Comparison for $(String(score_col))", :labels => labels, :linewidth => 1.5)
    merged_kwargs = merge(default_kwargs, kwargs)
    return plot_efdr_vs_qval(df, qval_col, efdr_cols; merged_kwargs...)
end

function plot_multiple_efdr_comparisons(df::DataFrame, score_qval_pairs::Vector{Tuple{Symbol,Symbol}}; method_types::Vector=[CombinedEFDR, PairedEFDR], layout=nothing, kwargs...)
    n_plots = length(score_qval_pairs)
    if isnothing(layout)
        n_cols = ceil(Int, sqrt(n_plots))
        n_rows = ceil(Int, n_plots / n_cols)
        layout = (n_rows, n_cols)
    end
    plots = []
    for (score_col, qval_col) in score_qval_pairs
        p = plot_efdr_comparison(df, score_col, qval_col; method_types=method_types, legend=(length(plots) == 0 ? :bottomright : false), kwargs...)
        push!(plots, p)
    end
    return plot(plots..., layout=layout, size=(800 * layout[2], 600 * layout[1]))
end

function save_efdr_plots(df::DataFrame, output_dir::String; score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)], method_types::Vector=[CombinedEFDR, PairedEFDR], formats::Vector{Symbol}=[:png, :pdf], title_suffix::AbstractString="")
    mkpath(output_dir)
    for (score_col, qval_col) in score_qval_pairs
        title = isempty(title_suffix) ? nothing : "EFDR Comparison for $(String(score_col)) $(title_suffix)"
        kwargs = isempty(title_suffix) ? (;) : (; title=title)
        p = plot_efdr_comparison(df, score_col, qval_col; method_types=method_types, kwargs...)
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_$(score_col).$(format)")
            savefig(p, filename)
            println("Saved: $filename")
        end
    end
    if length(score_qval_pairs) > 1
        p_combined = plot_multiple_efdr_comparisons(df, score_qval_pairs; method_types=method_types)
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_all.$(format)")
            savefig(p_combined, filename)
            println("Saved: $filename")
        end
    end
end

"""
Plot EFDR comparison across multiple replicate DataFrames on a single plot.
Colors are per method; linestyles per replicate.
"""
function plot_efdr_comparison_replicates(dfs::Vector{DataFrame}, score_col::Symbol, qval_col::Symbol;
                                         method_types::Vector=[CombinedEFDR, PairedEFDR],
                                         replicate_labels::Vector{String}=String[],
                                         method_colors=Dict(CombinedEFDR=>:blue, PairedEFDR=>:red),
                                         title="Entrapment vs Decoy FDR (Replicates)",
                                         legend=:bottomright,
                                         linewidth::Real=1.5)
    # Validate labels
    if !isempty(replicate_labels) && length(replicate_labels) != length(dfs)
        error("replicate_labels length must match dfs length")
    end
    # Determine auto-limits across all replicates/methods
    efdr_cols = Symbol[]
    for method_type in method_types
        method_name = method_type == CombinedEFDR ? "combined" : method_type == PairedEFDR ? "paired" : error("Unknown method")
        push!(efdr_cols, Symbol(String(score_col) * "_" * method_name * "_efdr"))
    end
    # Combine limits
    global_x = (0.0, 0.01)
    global_y = (0.0, 0.01)
    for df in dfs
        ax_x, ax_y = _auto_axis_limits(df, qval_col, efdr_cols)
        global_x = (0.0, max(global_x[2], ax_x[2]))
        global_y = (0.0, max(global_y[2], ax_y[2]))
    end
    p = plot(title=title, xlabel="Decoy FDR", ylabel="Entrapment FDR", xlims=global_x, ylims=global_y, legend=legend, size=(700, 550), dpi=300, fontfamily="Helvetica")
    max_val = min(global_x[2], global_y[2])
    plot!(p, [0, max_val], [0, max_val], label="y=x", linestyle=:dash, color=:gray, alpha=0.5)

    shown_label = Dict{UnionAll,Bool}(CombinedEFDR=>false, PairedEFDR=>false)
    for (rid, df) in enumerate(dfs)
        sorted_indices = sortperm(df[!, qval_col])
        sdf = df[sorted_indices, :]
        for method_type in method_types
            method_name = method_type == CombinedEFDR ? "combined" : "paired"
            efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
            if !hasproperty(sdf, efdr_col)
                @warn "Column $efdr_col not found in a replicate; skipping."
                continue
            end
            color = get(method_colors, method_type, :black)
            # Only one legend entry per method; all lines solid and slightly thinner
            label = shown_label[method_type] ? "" : titlecase(method_name)
            plot!(p, sdf[!, qval_col], sdf[!, efdr_col]; label=label, color=color, linestyle=:solid, linewidth=linewidth)
            shown_label[method_type] = true
        end
    end
    return p
end

function save_efdr_replicate_plots(dfs::Vector{DataFrame}, output_dir::String;
                                   score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                                   method_types::Vector=[CombinedEFDR, PairedEFDR],
                                   replicate_labels::Vector{String}=String[],
                                   formats::Vector{Symbol}=[:png, :pdf],
                                   title_suffix::AbstractString="")
    mkpath(output_dir)
    for (score_col, qval_col) in score_qval_pairs
        title = isempty(title_suffix) ? nothing : "Entrapment vs Decoy FDR (Replicates) $(title_suffix)"
        kwargs = isempty(title_suffix) ? (;) : (; title=title)
        p = plot_efdr_comparison_replicates(dfs, score_col, qval_col; method_types=method_types, replicate_labels=replicate_labels, kwargs...)
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_replicates_$(score_col).$(format)")
            savefig(p, filename)
            println("Saved: $filename")
        end
    end
end

function save_efdr_replicate_plots(pairdfs::Dict{Symbol, Vector{DataFrame}}, output_dir::String;
                                   score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                                   method_types::Vector=[CombinedEFDR, PairedEFDR],
                                   replicate_labels_map::Dict{Symbol, Vector{String}}=Dict{Symbol, Vector{String}}(),
                                   formats::Vector{Symbol}=[:png, :pdf],
                                   title_suffix::AbstractString="")
    mkpath(output_dir)
    for (score_col, qval_col) in score_qval_pairs
        dfs = get(pairdfs, score_col, DataFrame[])
        if isempty(dfs)
            @warn "No replicate data available for $(score_col); skipping."
            continue
        end
        labels = get(replicate_labels_map, score_col, String[])
        title = isempty(title_suffix) ? nothing : "Entrapment vs Decoy FDR (Replicates) $(title_suffix)"
        kwargs = isempty(title_suffix) ? (;) : (; title=title)
        p = plot_efdr_comparison_replicates(dfs, score_col, qval_col; method_types=method_types, replicate_labels=labels, kwargs...)
        for format in formats
            filename = joinpath(output_dir, "efdr_comparison_replicates_$(score_col).$(format)")
            savefig(p, filename)
            println("Saved: $filename")
        end
    end
end

export plot_efdr_vs_qval, plot_efdr_comparison, plot_multiple_efdr_comparisons, save_efdr_plots,
       plot_efdr_comparison_replicates, save_efdr_replicate_plots
