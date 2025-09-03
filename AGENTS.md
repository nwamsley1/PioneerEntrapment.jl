# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Julia source. Entry `PioneerEntrapment.jl` (exports/public API), `api.jl` (high-level workflows), `cli.jl` (CLI parser), `core/` (EFDR + pairing logic), `analysis/` (analysis helpers), `plotting/` (plots + save helpers).
- `test/`: `runtests.jl` with focused suites (pairing, protein EFDR). Add new tests as `test_*.jl` and include from `runtests.jl`.
- `bin/`: Executables `pioneer-entrapment` and `pioneer-entrapment-batch`.
- `docs/`: Documenter.jl setup (`docs/make.jl`, `docs/src/`).

## Build, Test, and Development Commands
- Install/instantiate: `julia --project -e 'using Pkg; Pkg.instantiate()'`
- Run tests: `julia --project -e 'using Pkg; Pkg.test()'`
- Coverage (local): `julia --project -e 'using Pkg; Pkg.test(; coverage=true)'`
- Build docs: `julia --project=docs docs/make.jl`
- CLI (single run): `bin/pioneer-entrapment --mode protein --protein-results path/to/proteins.arrow --outdir out`
- Batch: `bin/pioneer-entrapment-batch --base-dir <DIR> --library <LIB> --out-name efdr_output --dry-run`

## Coding Style & Naming Conventions
- Julia style: 4-space indent; no tabs or trailing whitespace.
- Names: Modules/Types `CamelCase`, functions/variables `lower_snake_case` (e.g., `add_protein_efdr_columns!`).
- Exports: update in `src/PioneerEntrapment.jl`; keep API surface tidy.
- Docstrings: triple-quoted with a brief example and argument notes.
- Formatting: prefer `JuliaFormatter.jl` (optional). Example: `julia -e 'using JuliaFormatter; format("src"); format("test")'`

## Testing Guidelines
- Framework: `Test`. Keep tests deterministic and small; avoid external files.
- Layout: add files like `test_feature_x.jl` and include from `runtests.jl`.
- Expectations: add tests for new functions and for bug fixes; maintain existing coverage.
- Run selectively in REPL: `] test PioneerEntrapment` or run files via `include("test/foo.jl")` inside a `@testset`.

## Commit & Pull Request Guidelines
- Commits: imperative, concise subject (<= 72 chars), scoped body when needed. Example: `fix(core): handle missing entrap_pair_id`.
- PRs: include description, rationale, and sample command/outputs; link issues; note breaking changes; update docs/exports and add tests.
- CI/docs: ensure `Pkg.test()` passes and `docs/make.jl` builds locally before requesting review.

## Security & Configuration Tips
- Large data: do not commit results; use paths outside the repo.
- Environment: scripts honor `JULIA_PROJECT` and `JULIA_DEPOT_PATH` (bin sets sensible defaults). Keep secrets out of flags and VCS.
