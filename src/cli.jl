"""
Simple CLI entrypoint. Provides a `julia_main(args)` that can be used by a bin script.

Usage examples:
  pioneer-entrapment --mode precursor \
    --precursor-results path/to/prec.arrow \
    --library path/to/library.arrow \
    --outdir out

  pioneer-entrapment --mode protein \
    --protein-results path/to/protein.arrow \
    --outdir out

  pioneer-entrapment --mode both \
    --precursor-results path/to/prec.arrow \
    --library path/to/library.arrow \
    --protein-results path/to/protein.arrow \
    --outdir out
"""

# Optional parsers for replicates config

function _parse_args(args)::Dict{String,Any}
    d = Dict{String,Any}()
    i = 1
    function take()
        i += 1
        if i > length(args)
            error("Missing value for flag")
        end
        return args[i]
    end
    while i <= length(args)
        a = args[i]
        if a == "--mode"
            d["mode"] = take()
        elseif a == "--precursor-results"
            d["precursor_results"] = take()
        elseif a == "--protein-results"
            d["protein_results"] = take()
        elseif a == "--library"
            d["library"] = take()
        elseif a == "--outdir"
            d["outdir"] = take()
        elseif a == "--r-lib"
            d["r_lib"] = parse(Float64, take())
        elseif a == "--plot-formats"
            # comma-separated list, e.g. png,pdf
            d["plot_formats"] = Symbol.(split(take(), ","))
        elseif a == "--paired-step" || a == "--paired-stride"
            d["paired_step"] = parse(Int, take())
        elseif a == "--replicates-config"
            d["replicates_config"] = take()
        elseif a == "--verbose"
            d["verbose"] = true
        elseif a == "-h" || a == "--help"
            println(usage())
            return Dict("__exit__" => 0)
        else
            error("Unknown argument: $a")
        end
        i += 1
    end
    return d
end

function usage()
    return """
PioneerEntrapment CLI

Required (depending on mode):
  --mode {precursor|protein|both}
  --precursor-results PATH         (for precursor or both)
  --library PATH                   (for precursor or both)
  --protein-results PATH           (for protein or both)
  --replicates-config PATH         (JSON/YAML/TOML list for replicates mode)

Optional:
  --outdir PATH                    (default: efdr_out)
  --r-lib FLOAT                    (default: 1.0)
  --paired-step INT                (default: 5) stride for paired EFDR
  --plot-formats LIST              (e.g., png,pdf)
  --verbose                        (enable verbose logging)
  -h, --help
"""
end

# Load replicate spec file: supports JSON (.json), YAML (.yml/.yaml), or TOML (.toml)
function _load_replicates_config(path::AbstractString)
    lower = lowercase(path)
    if endswith(lower, ".json")
        if isdefined(Main, :JSON)
            data = Main.JSON.parsefile(path)
        else
            error("JSON.jl not available. Install with: using Pkg; Pkg.add(\"JSON\"); then `using JSON` before calling, or use a TOML config.")
        end
    elseif endswith(lower, ".yaml") || endswith(lower, ".yml")
        if isdefined(Main, :YAML)
            data = Main.YAML.load_file(path)
        else
            error("YAML.jl not available. Install with: using Pkg; Pkg.add(\"YAML\"); then `using YAML` before calling, or use a TOML config.")
        end
    elseif endswith(lower, ".toml")
        data = Base.TOML.parsefile(path)
        data = get(data, "replicates", data)
    else
        error("Unsupported config extension for $path. Use .json, .yaml/.yml, or .toml")
    end
    # Normalize to vector of dict-like entries
    reps = Vector{Any}(data)
    normalized = Vector{NamedTuple}(undef, length(reps))
    for (i, r) in enumerate(reps)
        # Allow keys with different casings
        pr = get(r, "precursor_results_path", get(r, :precursor_results_path, nothing))
        lib = get(r, "library_precursors_path", get(r, :library_precursors_path, nothing))
        lbl = get(r, "label", get(r, :label, "rep$(i)"))
        if pr === nothing || lib === nothing
            error("Replicate #$i is missing required precursor_results_path or library_precursors_path")
        end
        normalized[i] = (precursor_results_path=String(pr), library_precursors_path=String(lib), label=String(lbl))
    end
    return normalized
end

function julia_main(args)::Int
    try
        parsed = _parse_args(args)
        if haskey(parsed, "__exit__")
            return parsed["__exit__"]
        end
        mode = get(parsed, "mode", "")
        outdir = get(parsed, "outdir", "efdr_out")
        r_lib = get(parsed, "r_lib", 1.0)
        plot_formats = get(parsed, "plot_formats", Symbol[:png, :pdf])
        paired_step = get(parsed, "paired_step", 5)
        verbose = get(parsed, "verbose", true)

        if mode == "replicates" || haskey(parsed, "replicates_config")
            if !haskey(parsed, "replicates_config")
                println("Missing required --replicates-config for replicates mode\n\n" * usage())
                return 2
            end
            reps = _load_replicates_config(parsed["replicates_config"])
            run_efdr_replicate_plots(reps; output_dir=outdir, r_lib=r_lib, paired_stride=paired_step, plot_formats=plot_formats, verbose=verbose)
            return 0
        end

        if mode == "precursor"
            pr = get(parsed, "precursor_results", nothing)
            lib = get(parsed, "library", nothing)
            if pr === nothing || lib === nothing
                println("Missing required --precursor-results or --library\n\n" * usage())
                return 2
            end
            run_efdr_analysis(pr, lib; output_dir=outdir, r_lib=r_lib, paired_stride=paired_step, plot_formats=plot_formats, verbose=verbose)
            return 0
        elseif mode == "protein"
            prot = get(parsed, "protein_results", nothing)
            if prot === nothing
                println("Missing required --protein-results\n\n" * usage())
                return 2
            end
            run_protein_efdr_analysis(prot; output_dir=outdir, r_lib=r_lib, paired_stride=paired_step, plot_formats=plot_formats, verbose=verbose)
            return 0
        elseif mode == "both"
            pr = get(parsed, "precursor_results", nothing)
            lib = get(parsed, "library", nothing)
            prot = get(parsed, "protein_results", nothing)
            if pr === nothing || lib === nothing || prot === nothing
                println("Missing required inputs for both mode\n\n" * usage())
                return 2
            end
            run_both_analyses(; precursor_results_path=pr, library_precursors_path=lib, protein_results_path=prot, output_dir=outdir, r_lib=r_lib, paired_stride=paired_step, plot_formats=plot_formats, verbose=verbose)
            return 0
        else
            println("Invalid or missing --mode\n\n" * usage())
            return 2
        end
    catch err
        println("Error: ", err)
        println(usage())
        return 1
    end
end

export julia_main
