using Arrow
using DataFrames
using Printf
using Markdown

"""
Load a tabular file into a DataFrame.
Supports Arrow/Feather (`.arrow`, `.feather`) and delimited text (`.csv`, `.tsv`, `.txt`).
Also accepts a directory path (e.g., a `.poin` folder) and will locate a suitable table inside,
preferring `precursors_table.arrow`, then any `.arrow`/`.feather`, then `.tsv`/`.csv`.
Falls back to trying Arrow, then CSV regardless of extension.
"""
function _load_table(path::AbstractString)
    # If given a directory (e.g., a .poin bundle), resolve to a concrete file inside
    if isdir(path)
        entries = readdir(path)
        # Priority list
        preferred = [
            "precursors_table.arrow",
            "library_precursors.arrow",
            "precursors.arrow",
        ]
        for name in preferred
            if name in entries
                path = joinpath(path, name)
                break
            end
        end
        # If no preferred name matched, pick first Arrow/Feather, else TSV/CSV
        if isdir(path)
            # still a dir, so choose by extension
            arrow_like = filter(x -> endswith(lowercase(x), ".arrow") || endswith(lowercase(x), ".feather"), entries)
            if !isempty(arrow_like)
                path = joinpath(path, first(arrow_like))
            else
                txt_like = filter(x -> endswith(lowercase(x), ".tsv") || endswith(lowercase(x), ".csv"), entries)
                if !isempty(txt_like)
                    path = joinpath(path, first(txt_like))
                end
            end
        end
    end

    lower = lowercase(path)
    # Prefer extension-based dispatch where possible
    if endswith(lower, ".arrow") || endswith(lower, ".feather")
        try
            return DataFrame(Arrow.Table(path); copycols=true)
        catch
            try
                open(path, "r") do io
                    return DataFrame(Arrow.Table(io; mmap=false); copycols=true)
                end
            catch
                return DataFrame(Arrow.Table(path; mmap=false); copycols=true)
            end
        end
    elseif endswith(lower, ".csv") || endswith(lower, ".tsv") || endswith(lower, ".txt")
        delim = endswith(lower, ".tsv") ? '\t' : ','
        return DataFrame(CSV.File(path; delim=delim))
    else
        # try Arrow, then CSV
        try
            return DataFrame(Arrow.Table(path); copycols=true)
        catch
            try
                return DataFrame(Arrow.Table(path; mmap=false); copycols=true)
            catch
                return DataFrame(CSV.File(path))
            end
        end
    end
end

"""
    run_efdr_plots(results_dir::String, library_path::String; output_dir=joinpath(results_dir, "efdr_out"), r_lib=1.0, paired_stride=5, plot_formats=[:png,:pdf], verbose=true)

Convenience entry point that looks for standard filenames in `results_dir`:
- precursors_long.arrow (or .tsv) for precursor-level
- protein_groups_long.arrow (or .tsv) for protein-level
Runs the appropriate analyses, writing outputs into `output_dir` (or subfolders if both).
"""
function run_efdr_plots(results_dir::String, library_path::String;
                        output_dir::String=joinpath(results_dir, "efdr_out"),
                        r_lib::Float64=1.0,
                        paired_stride::Int=5,
                        plot_formats::Vector{Symbol} = [:png, :pdf],
                        verbose::Bool=true)
    prec = if isfile(joinpath(results_dir, "precursors_long.arrow"))
        joinpath(results_dir, "precursors_long.arrow")
    elseif isfile(joinpath(results_dir, "precursors_long.tsv"))
        joinpath(results_dir, "precursors_long.tsv")
    else
        nothing
    end
    prot = if isfile(joinpath(results_dir, "protein_groups_long.arrow"))
        joinpath(results_dir, "protein_groups_long.arrow")
    elseif isfile(joinpath(results_dir, "protein_groups_long.tsv"))
        joinpath(results_dir, "protein_groups_long.tsv")
    else
        nothing
    end
    if prec !== nothing && prot !== nothing
        return run_both_analyses(; precursor_results_path=prec,
                                    library_precursors_path=library_path,
                                    protein_results_path=prot,
                                    output_dir=output_dir,
                                    r_lib=r_lib,
                                    paired_stride=paired_stride,
                                    plot_formats=plot_formats,
                                    verbose=verbose)
    elseif prec !== nothing
        return run_efdr_analysis(prec, library_path; output_dir=output_dir, r_lib=r_lib, paired_stride=paired_stride, plot_formats=plot_formats, verbose=verbose)
    elseif prot !== nothing
        return run_protein_efdr_analysis(prot; output_dir=output_dir, r_lib=r_lib, paired_stride=paired_stride, plot_formats=plot_formats, verbose=verbose)
    else
        error("No standard result files found in $(results_dir). Expected precursors_long.(arrow|tsv) and/or protein_groups_long.(arrow|tsv)")
    end
end

"""
    run_efdr_replicate_plots(replicates; output_dir="efdr_out", score_qval_pairs=[(:global_prob, :global_qval), (:prec_prob, :qval)], r_lib=1.0, paired_stride=5, plot_formats=[:png,:pdf], verbose=true)

Compute EFDR for multiple (precursor_results_path, library_precursors_path) replicates and plot on shared figures.

`replicates` is a Vector of NamedTuples with fields:
- `precursor_results_path::String`
- `library_precursors_path::String`
- `label::String` (optional, used in plot legend)
"""
function run_efdr_replicate_plots(replicates::Vector; output_dir::String="efdr_out",
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                                  protein_score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
                                  r_lib::Float64=1.0, paired_stride::Int=5,
                                  plot_formats::Vector{Symbol}=[:png, :pdf], verbose::Bool=true)
    # Prepare per-score collections of replicate dataframes and labels
    pairdfs = Dict{Symbol, Vector{DataFrame}}()
    labels_map = Dict{Symbol, Vector{String}}()
    for (score_col, _) in score_qval_pairs
        pairdfs[score_col] = DataFrame[]
        labels_map[score_col] = String[]
    end
    for (score_col, _) in protein_score_qval_pairs
        pairdfs[score_col] = get(pairdfs, score_col, DataFrame[])
        labels_map[score_col] = get(labels_map, score_col, String[])
    end

    # Small helpers to work with NamedTuple or Dict-like entries
    hasprop(x, s::Symbol) = s in propertynames(x)
    getprop(x, s::Symbol, default) = hasprop(x, s) ? getfield(x, s) : default

    for (idx, rep) in enumerate(replicates)
        # Support NamedTuple entries with fields: :precursor_results_path or :rep_dir, :library_precursors_path, optional :label
        pr_path = getprop(rep, :precursor_results_path, nothing)
        lib_path = getprop(rep, :library_precursors_path, nothing)
        label = getprop(rep, :label, "rep$(idx)")
        # Ensure rep_dir is defined before any use
        rep_dir = getprop(rep, :rep_dir, nothing)
        if rep_dir === nothing && pr_path !== nothing
            rep_dir = dirname(String(pr_path))
        end

        # Resolve paths from rep_dir if precursor path was not provided
        if pr_path === nothing && rep_dir !== nothing
            if isfile(joinpath(rep_dir, "precursors_long.arrow"))
                pr_path = joinpath(rep_dir, "precursors_long.arrow")
            elseif isfile(joinpath(rep_dir, "precursors_long.tsv"))
                pr_path = joinpath(rep_dir, "precursors_long.tsv")
            end
        end
        if pr_path === nothing
            error("Replicate #$(idx) missing precursor results: provide :precursor_results_path or ensure precursors_long.(arrow|tsv) exists in :rep_dir")
        end
        if lib_path === nothing
            error("Replicate #$(idx) missing :library_precursors_path")
        end

        verbose && println("[rep$(idx)] Loading data...")
        prec_results = _load_table(String(pr_path))
        library_precursors = _load_table(String(lib_path))

        # Filter non-targets if present
        if hasproperty(prec_results, :target)
            filter!(x -> x.target, prec_results)
        end

        if !hasproperty(library_precursors, :mod_key)
            library_precursors[!, :mod_key] = map(x -> getModKey(x), library_precursors.structural_mods)
        end
        assign_entrapment_pairs!(library_precursors)
        add_entrap_pair_ids!(prec_results, library_precursors)

        # Identify which requested pairs are "global" vs per-file (precursor)
        global_pairs = [(s,q) for (s,q) in score_qval_pairs if occursin("global", String(s))]
        perfile_pairs = [(s,q) for (s,q) in score_qval_pairs if !occursin("global", String(s))]

        # Per-file EFDRs
        if !isempty(perfile_pairs)
            add_original_target_scores!(prec_results, library_precursors, [s for (s,_) in perfile_pairs])
            add_efdr_columns!(prec_results, library_precursors; score_qval_pairs=perfile_pairs, r=r_lib, paired_stride=paired_stride)
            for (s, _) in perfile_pairs
                # Push only if EFDR columns exist
                for method_type in (CombinedEFDR, PairedEFDR)
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(s) * "_" * method_name * "_efdr")
                    if hasproperty(prec_results, efdr_col)
                        push!(pairdfs[s], prec_results)
                        push!(labels_map[s], label)
                        break
                    end
                end
            end
        end

        # Global EFDRs (precursor)
        if !isempty(global_pairs)
            global_df = create_global_results_df(prec_results; score_col=first(global_pairs)[1])
            add_original_target_scores!(global_df, library_precursors, [s for (s,_) in global_pairs])
            add_efdr_columns!(global_df, library_precursors; score_qval_pairs=global_pairs, r=r_lib, paired_stride=paired_stride)
            for (s, _) in global_pairs
                for method_type in (CombinedEFDR, PairedEFDR)
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(s) * "_" * method_name * "_efdr")
                    if hasproperty(global_df, efdr_col)
                        push!(pairdfs[s], global_df)
                        push!(labels_map[s], label)
                        break
                    end
                end
            end
        end

        # Protein-level EFDRs (if protein file exists alongside the replicate)
        # Infer a protein file path from the replicate entry if possible
        # rep_dir already established above
        if rep_dir !== nothing
            prot_path = if isfile(joinpath(rep_dir, "protein_groups_long.arrow"))
                joinpath(rep_dir, "protein_groups_long.arrow")
            elseif isfile(joinpath(rep_dir, "protein_groups_long.tsv"))
                joinpath(rep_dir, "protein_groups_long.tsv")
            else
                nothing
            end
            if prot_path !== nothing
                # Load and compute protein EFDR columns
                protein_results = _load_table(prot_path)
                # Separate global/per-file
                prot_global_pairs = [(s,q) for (s,q) in protein_score_qval_pairs if occursin("global", String(s))]
                prot_perfile_pairs = [(s,q) for (s,q) in protein_score_qval_pairs if !occursin("global", String(s))]
                if !isempty(prot_perfile_pairs)
                    add_original_target_protein_scores!(protein_results, [s for (s,_) in prot_perfile_pairs])
                    add_protein_efdr_columns!(protein_results; score_qval_pairs=prot_perfile_pairs, r=r_lib, paired_stride=paired_stride)
                    for (s, _) in prot_perfile_pairs
                        for method_type in (CombinedEFDR, PairedEFDR)
                            method_name = method_type == CombinedEFDR ? "combined" : "paired"
                            efdr_col = Symbol(String(s) * "_" * method_name * "_efdr")
                            if hasproperty(protein_results, efdr_col)
                                push!(pairdfs[s], protein_results)
                                push!(labels_map[s], label)
                                break
                            end
                        end
                    end
                end
                if !isempty(prot_global_pairs)
                    global_prot_df = create_global_protein_results_df(protein_results; score_col=first(prot_global_pairs)[1])
                    add_original_target_protein_scores!(global_prot_df, [s for (s,_) in prot_global_pairs])
                    add_protein_efdr_columns!(global_prot_df; score_qval_pairs=prot_global_pairs, r=r_lib, paired_stride=paired_stride)
                    for (s, _) in prot_global_pairs
                        for method_type in (CombinedEFDR, PairedEFDR)
                            method_name = method_type == CombinedEFDR ? "combined" : "paired"
                            efdr_col = Symbol(String(s) * "_" * method_name * "_efdr")
                            if hasproperty(global_prot_df, efdr_col)
                                push!(pairdfs[s], global_prot_df)
                                push!(labels_map[s], label)
                                break
                            end
                        end
                    end
                end
            end
        end
    end

    # Save per-pair replicate plots
    save_efdr_replicate_plots(pairdfs, output_dir; score_qval_pairs=score_qval_pairs, replicate_labels_map=labels_map, formats=plot_formats)
    # Also save plots for any protein score pairs that have data
    save_efdr_replicate_plots(pairdfs, output_dir; score_qval_pairs=protein_score_qval_pairs, replicate_labels_map=labels_map, formats=plot_formats)
    return (
        dataframes_by_score = pairdfs,
        labels_by_score = labels_map,
        output_dir = output_dir,
    )
end

"""
    run_efdr_analysis(prec_results_path::String, library_precursors_path::String;
                      output_dir::String="efdr_out",
                      method_types=[CombinedEFDR, PairedEFDR],
                      score_qval_pairs=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                      r_lib::Float64=1.0,
                      plot_formats=[:png, :pdf],
                      verbose::Bool=true)

Run empirical FDR analysis on precursor-level data with entrapment sequences.
Accepts Arrow/CSV inputs for convenience.
"""
function run_efdr_analysis(prec_results_path::String, library_precursors_path::String;
                          output_dir::String="efdr_out",
                          method_types::Vector=[CombinedEFDR, PairedEFDR],
                          score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_prob, :global_qval), (:prec_prob, :qval)],
                          r_lib::Float64=1.0,
                          paired_stride::Int=5,
                          plot_formats::Vector{Symbol}=[:png, :pdf],
                          verbose::Bool=true)

    mkpath(output_dir)
    output_files = String[]

    verbose && println("Loading data...")
    prec_results = _load_table(prec_results_path)
    library_precursors = _load_table(library_precursors_path)

    verbose && println("Loaded $(nrow(prec_results)) precursor results")
    verbose && println("Loaded $(nrow(library_precursors)) library precursors")

    original_rows = nrow(prec_results)
    if hasproperty(prec_results, :target)
        non_target_rows = sum(.!prec_results.target)
        if non_target_rows > 0
            verbose && @warn "Filtering out $non_target_rows non-target (decoy) rows before EFDR calculation"
            filter!(x -> x.target, prec_results)
        end
    else
        verbose && @warn "No 'target' column found. Assuming all rows are targets."
    end

    if !hasproperty(library_precursors, :mod_key)
        verbose && println("Adding mod_key column to library...")
        library_precursors[!, :mod_key] = map(x -> getModKey(x), library_precursors.structural_mods)
    end

    verbose && println("Assigning entrapment pairs...")
    assign_entrapment_pairs!(library_precursors)

    add_entrap_pair_ids!(prec_results, library_precursors)

    verbose && println("Separating global and per-file analyses...")
    global_scores = Symbol[]
    perfile_scores = Symbol[]
    for (score_col, _) in score_qval_pairs
        if occursin("global", String(score_col))
            push!(global_scores, score_col)
        else
            push!(perfile_scores, score_col)
        end
    end

    global_results_df = nothing
    if !isempty(global_scores)
        verbose && println("Creating global results dataframe...")
        global_results_df = create_global_results_df(prec_results; score_col=global_scores[1])
        verbose && println("Global dataframe has $(nrow(global_results_df)) unique precursors")
    end

    if !isempty(perfile_scores)
        verbose && println("Processing per-file scores...")
        add_original_target_scores!(prec_results, library_precursors, perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        if !isempty(perfile_pairs)
            add_efdr_columns!(prec_results, library_precursors;
                              method_types=method_types,
                              score_qval_pairs=perfile_pairs,
                              r=r_lib,
                              paired_stride=paired_stride)
        end
    end

    if !isnothing(global_results_df) && !isempty(global_scores)
        verbose && println("Processing global scores...")
        add_original_target_scores!(global_results_df, library_precursors, global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        if !isempty(global_pairs)
            add_efdr_columns!(global_results_df, library_precursors;
                              method_types=method_types,
                              score_qval_pairs=global_pairs,
                              r=r_lib,
                              paired_stride=paired_stride)
        end
    end

    comparison_results = Dict{Tuple{Symbol,Symbol}, DataFrame}()
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            verbose && println("Comparing EFDR methods for $score_col/$qval_col (per-file)...")
            comparison_df = compare_efdr_methods(prec_results, qval_col, score_col, library_precursors)
            comparison_results[(score_col, qval_col)] = comparison_df
        end
    end
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                verbose && println("Comparing EFDR methods for $score_col/$qval_col (global)...")
                comparison_df = compare_efdr_methods(global_results_df, qval_col, score_col, library_precursors)
                comparison_results[(score_col, qval_col)] = comparison_df
            end
        end
    end

    calibration_results = Dict{Symbol, Tuple{DataFrame, Float64}}()
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            for method_type in method_types
                method_name = method_type == CombinedEFDR ? "combined" : "paired"
                efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                verbose && println("Calculating calibration error for $efdr_col (per-file)...")
                cal_data, cal_error = calculate_efdr_calibration_error(
                    prec_results, qval_col, efdr_col, library_precursors
                )
                calibration_results[efdr_col] = (cal_data, cal_error)
            end
        end
    end
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                for method_type in method_types
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                    verbose && println("Calculating calibration error for $efdr_col (global)...")
                    cal_data, cal_error = calculate_efdr_calibration_error(
                        global_results_df, qval_col, efdr_col, library_precursors
                    )
                    calibration_results[efdr_col] = (cal_data, cal_error)
                end
            end
        end
    end

    verbose && println("Generating plots...")
    if !isempty(perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        save_efdr_plots(prec_results, output_dir; score_qval_pairs=perfile_pairs, method_types=method_types, formats=plot_formats)
    end
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        save_efdr_plots(global_results_df, output_dir; score_qval_pairs=global_pairs, method_types=method_types, formats=plot_formats)
    end
    for (score_col, _) in score_qval_pairs
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_$(score_col).$(format)"))
        end
    end
    if length(score_qval_pairs) > 1
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_all.$(format)"))
        end
    end

    verbose && println("Generating markdown report...")
    report_path = joinpath(output_dir, "efdr_analysis_report.md")
    report_df = if !isnothing(global_results_df) && isempty(perfile_scores)
        global_results_df
    else
        prec_results
    end
    generate_markdown_report(report_df, library_precursors, comparison_results,
                             calibration_results, score_qval_pairs, method_types,
                             original_rows, output_dir, report_path)
    push!(output_files, report_path)

    verbose && println("Saving processed data...")
    saved_dfs = Dict{String, DataFrame}()
    if !isempty(perfile_scores)
        data_path = joinpath(output_dir, "prec_results_with_efdr.arrow")
        Arrow.write(data_path, prec_results)
        push!(output_files, data_path)
        saved_dfs["per_file"] = prec_results
    end
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_data_path = joinpath(output_dir, "global_results_with_efdr.arrow")
        Arrow.write(global_data_path, global_results_df)
        push!(output_files, global_data_path)
        saved_dfs["global"] = global_results_df
    end

    verbose && println("\nAnalysis complete! Output saved to: $output_dir")

    return (
        filtered_data = length(saved_dfs) == 2 ? saved_dfs : get(saved_dfs, "global", get(saved_dfs, "per_file", DataFrame())),
        comparison_results = comparison_results,
        calibration_results = calibration_results,
        output_files = output_files
    )
end

"""
    generate_markdown_report(prec_results, library_precursors, comparison_results,
                           calibration_results, score_qval_pairs, method_types,
                           original_rows, output_dir, output_path)
"""
function generate_markdown_report(prec_results::DataFrame, library_precursors::DataFrame,
                                  comparison_results::Dict, calibration_results::Dict,
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}},
                                  method_types::Vector, original_rows::Int,
                                  output_dir::String, output_path::String)

    open(output_path, "w") do io
        println(io, "# Empirical FDR Analysis Report")
        println(io, "\nGenerated: $(now())")
        println(io)
        println(io, "## Data Summary")
        println(io)
        println(io, "- Original precursor results: $(original_rows)")
        filtered_rows = original_rows - nrow(prec_results)
        if filtered_rows > 0
            println(io, "- Filtered non-target rows: $(filtered_rows)")
        end
        println(io, "- Analyzed precursor results: $(nrow(prec_results))")
        println(io, "- Library precursors: $(nrow(library_precursors))")
        println(io)
        println(io, "## EFDR Comparisons")
        for ((score_col, qval_col), df) in comparison_results
            println(io, "### $(score_col) / $(qval_col)")
            println(io, "- Rows: $(nrow(df))")
        end
        println(io)
        println(io, "## Calibration")
        for (efdr_col, (_, cal_error)) in calibration_results
            @printf(io, "- %s: mean abs error = %.4f\n", String(efdr_col), cal_error)
        end
        println(io)
        println(io, "## Plots")
        for (score_col, _) in score_qval_pairs
            println(io, "- efdr_comparison_$(score_col).png")
        end
        if length(score_qval_pairs) > 1
            println(io, "- efdr_comparison_all.png")
        end
        println(io)
        println(io, "## Parameters")
        println(io, "- EFDR Methods: ", join([m == CombinedEFDR ? "Combined" : "Paired" for m in method_types], ", "))
        println(io, "- Score/Q-value pairs: ", join(["$s/$q" for (s,q) in score_qval_pairs], ", "))
        println(io, "- Output dir: `$output_dir`")
    end
end

"""
    run_protein_efdr_analysis(protein_results_path::String;
                              output_dir::String="efdr_out",
                              method_types=[CombinedEFDR, PairedEFDR],
                              score_qval_pairs=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
                              r_lib::Float64=1.0,
                              plot_formats=[:png, :pdf],
                              verbose::Bool=true)
"""
function run_protein_efdr_analysis(protein_results_path::String;
                                  output_dir::String="efdr_out",
                                  method_types::Vector=[CombinedEFDR, PairedEFDR],
                                  score_qval_pairs::Vector{Tuple{Symbol,Symbol}}=[(:global_pg_score, :global_qval), (:pg_score, :qval)],
                                  r_lib::Float64=1.0,
                                  paired_stride::Int=5,
                                  plot_formats::AbstractVector=[:png, :pdf],
                                  verbose::Bool=true)

    mkpath(output_dir)
    output_files = String[]

    verbose && println("Loading protein data...")
    protein_results = _load_table(protein_results_path)
    is_sorted = issorted(protein_results, [:pg_score, :entrap_id], rev = [true, false])
    @info is_sorted ? "Protein results are sorted by pg_score and entrap_id" :
        "Protein results are NOT sorted by pg_score and entrap_id, sorting now"
    sort!(protein_results, [:pg_score, :entrap_id], rev = [true, false])

    verbose && println("Loaded $(nrow(protein_results)) protein results")

    protein_cols = [:protein, :entrap_id, :pg_score]
    if !all(col -> hasproperty(protein_results, col), protein_cols)
        error("Input file does not appear to be protein data. Expected columns: $protein_cols")
    end

    original_rows = nrow(protein_results)
    if hasproperty(protein_results, :target)
        non_target_rows = sum(.!protein_results.target)
        if non_target_rows > 0
            verbose && @warn "Filtering out $non_target_rows non-target (decoy) rows before EFDR calculation"
            filter!(x -> x.target, protein_results)
        end
    else
        throw("No 'target' column found.")
    end

    verbose && println("Separating global and per-file analyses...")
    global_scores = Symbol[]
    perfile_scores = Symbol[]
    for (score_col, _) in score_qval_pairs
        if occursin("global", String(score_col))
            push!(global_scores, score_col)
        else
            push!(perfile_scores, score_col)
        end
    end

    global_results_df = nothing
    if !isempty(global_scores)
        verbose && println("Creating global protein results dataframe...")
        global_results_df = create_global_protein_results_df(protein_results; score_col=global_scores[1])
        verbose && println("Global dataframe has $(nrow(global_results_df)) unique proteins")
    end

    # Normalize formats to Vector{Symbol}
    plot_formats = isempty(plot_formats) ? Symbol[] : Symbol.(plot_formats)

    if !isempty(perfile_scores)
        verbose && println("Processing per-file protein scores...")
        add_original_target_protein_scores!(protein_results, perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        if !isempty(perfile_pairs)
            add_protein_efdr_columns!(protein_results; method_types=method_types, score_qval_pairs=perfile_pairs, r=r_lib, paired_stride=paired_stride)
        end
    end

    if !isnothing(global_results_df) && !isempty(global_scores)
        verbose && println("Processing global protein scores...")
        add_original_target_protein_scores!(global_results_df, global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        if !isempty(global_pairs)
            add_protein_efdr_columns!(global_results_df; method_types=method_types, score_qval_pairs=global_pairs, r=r_lib, paired_stride=paired_stride)
        end
    end

    comparison_results = Dict{Tuple{Symbol,Symbol}, DataFrame}()
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            verbose && println("Comparing EFDR methods for $score_col/$qval_col (per-file)...")
            comparison_df = compare_protein_efdr_methods(protein_results, qval_col, score_col)
            comparison_results[(score_col, qval_col)] = comparison_df
        end
    end
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                verbose && println("Comparing EFDR methods for $score_col/$qval_col (global)...")
                comparison_df = compare_protein_efdr_methods(global_results_df, qval_col, score_col)
                comparison_results[(score_col, qval_col)] = comparison_df
            end
        end
    end

    calibration_results = Dict{Symbol, Tuple{DataFrame, Float64}}()
    for (score_col, qval_col) in score_qval_pairs
        if score_col in perfile_scores
            for method_type in method_types
                method_name = method_type == CombinedEFDR ? "combined" : "paired"
                efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                if hasproperty(protein_results, efdr_col)
                    verbose && println("Calculating calibration error for $efdr_col (per-file)...")
                    cal_data, cal_error = calculate_protein_efdr_calibration_error(
                        protein_results, qval_col, efdr_col
                    )
                    calibration_results[efdr_col] = (cal_data, cal_error)
                end
            end
        end
    end
    if !isnothing(global_results_df)
        for (score_col, qval_col) in score_qval_pairs
            if score_col in global_scores
                for method_type in method_types
                    method_name = method_type == CombinedEFDR ? "combined" : "paired"
                    efdr_col = Symbol(String(score_col) * "_" * method_name * "_efdr")
                    if hasproperty(global_results_df, efdr_col)
                        verbose && println("Calculating calibration error for $efdr_col (global)...")
                        cal_data, cal_error = calculate_protein_efdr_calibration_error(
                            global_results_df, qval_col, efdr_col
                        )
                        calibration_results[efdr_col] = (cal_data, cal_error)
                    end
                end
            end
        end
    end

    verbose && println("Generating plots...")
    if !isempty(perfile_scores)
        perfile_pairs = [(s, q) for (s, q) in score_qval_pairs if s in perfile_scores]
        save_efdr_plots(protein_results, output_dir; score_qval_pairs=perfile_pairs, method_types=method_types, formats=plot_formats)
    end
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_pairs = [(s, q) for (s, q) in score_qval_pairs if s in global_scores]
        save_efdr_plots(global_results_df, output_dir; score_qval_pairs=global_pairs, method_types=method_types, formats=plot_formats)
    end
    for (score_col, _) in score_qval_pairs
        for format in plot_formats
            push!(output_files, joinpath(output_dir, "efdr_comparison_$(score_col).$(format)"))
        end
    end

    verbose && println("Generating markdown report...")
    report_path = joinpath(output_dir, "protein_efdr_analysis_report.md")
    open(report_path, "w") do io
        println(io, "# Protein-Level Empirical FDR Analysis Report")
        println(io, "\nGenerated: $(now())")
        println(io, "\n## Data Summary")
        println(io, "- Original protein results: $(original_rows)")
        println(io, "- Analyzed protein results: $(nrow(protein_results))")
    end
    push!(output_files, report_path)

    verbose && println("Saving processed data...")
    saved_dfs = Dict{String, DataFrame}()
    if !isempty(perfile_scores)
        data_path = joinpath(output_dir, "protein_results_with_efdr.arrow")
        Arrow.write(data_path, protein_results)
        push!(output_files, data_path)
        saved_dfs["per_file"] = protein_results
    end
    if !isnothing(global_results_df) && !isempty(global_scores)
        global_data_path = joinpath(output_dir, "global_protein_results_with_efdr.arrow")
        Arrow.write(global_data_path, global_results_df)
        push!(output_files, global_data_path)
        saved_dfs["global"] = global_results_df
    end

    verbose && println("\nProtein analysis complete! Output saved to: $output_dir")

    return (
        filtered_data = length(saved_dfs) == 2 ? saved_dfs : get(saved_dfs, "global", get(saved_dfs, "per_file", DataFrame())),
        comparison_results = comparison_results,
        calibration_results = calibration_results,
        output_files = output_files
    )
end

"""
    run_both_analyses(; precursor_results_path::AbstractString,
                        library_precursors_path::AbstractString,
                        protein_results_path::AbstractString,
                        output_dir::AbstractString = "efdr_out",
                        r_lib::Float64 = 1.0,
                        plot_formats::Vector{Symbol} = [:png, :pdf],
                        verbose::Bool = true)

Run both the precursor-level and protein-level analyses. Returns a NamedTuple
with both results, writing outputs into `output_dir/precursor` and `output_dir/protein`.
"""
function run_both_analyses(; precursor_results_path::AbstractString,
                        library_precursors_path::AbstractString,
                        protein_results_path::AbstractString,
                        output_dir::AbstractString = "efdr_out",
                        r_lib::Float64 = 1.0,
                        paired_stride::Int = 5,
                        plot_formats::Vector{Symbol} = [:png, :pdf],
                        verbose::Bool = true)

    out_prec = joinpath(output_dir, "precursor")
    out_prot = joinpath(output_dir, "protein")

    prec = run_efdr_analysis(precursor_results_path, library_precursors_path;
                             output_dir=out_prec, r_lib=r_lib, paired_stride=paired_stride, plot_formats=plot_formats, verbose=verbose)
    prot = run_protein_efdr_analysis(protein_results_path;
                                     output_dir=out_prot, r_lib=r_lib, paired_stride=paired_stride, plot_formats=plot_formats, verbose=verbose)

    return (precursor=prec, protein=prot)
end
