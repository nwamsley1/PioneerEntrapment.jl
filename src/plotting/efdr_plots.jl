using Plots
using DataFrames
using Statistics

function plot_efdr_vs_qval(df::DataFrame, qval_col::Symbol, efdr_cols::Vector{Symbol};
                          title="Empirical FDR vs Q-value",
                          xlabel="Q-value",
                          ylabel="Empirical FDR",
                          labels=nothing,
                          colors=nothing,
                          xlims=(0, 0.015),
                          ylims=(0, 0.015),
                          legend=:bottomright,
                          diagonal=true)
    sorted_indices = sortperm(df[!, qval_col])
    sorted_df = df[sorted_indices, :]
    p = plot(title=title, xlabel=xlabel, ylabel=ylabel, xlims=xlims, ylims=ylims, legend=legend, size=(600, 500), dpi=300)
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
        plot!(p, sorted_df[!, qval_col], sorted_df[!, efdr_col], label=label, color=color, linewidth=2, alpha=0.8)
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
    default_kwargs = Dict(:title => "EFDR Comparison for $(String(score_col))", :labels => labels)
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

function save_efdr_plots(df::DataFrame, output_dir::String; score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)], method_types::Vector=[CombinedEFDR, PairedEFDR], formats::Vector{Symbol}=[:png, :pdf])
    mkpath(output_dir)
    for (score_col, qval_col) in score_qval_pairs
        p = plot_efdr_comparison(df, score_col, qval_col; method_types=method_types)
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

export plot_efdr_vs_qval, plot_efdr_comparison, plot_multiple_efdr_comparisons, save_efdr_plots

