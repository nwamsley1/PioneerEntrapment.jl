# Entrapment Pairing and EFDR Computation

This page explains how PioneerEntrapment.jl defines target/entrapment groupings, assigns pair identifiers, derives “original target” scores, and computes empirical FDR (EFDR) using both the Combined and Paired approaches.

## Key Concepts

- Target vs Entrapment: Original targets are labeled with `0` and entrapments with non-zero integers.
  - Precursor level: `:entrapment_group_id` (from the spectral library).
  - Protein level: `:entrap_id` (from the protein results/library).
- Entrapment Pair: A stable identifier that groups an original target with its corresponding entrapment analog(s). Pairs allow paired EFDR to reason about target–entrapment relationships within the same biological entity.
- Original Target Score: For a given entrapment row, the score of its corresponding original target (its pair-mate) at the same file (and species, for proteins). This is materialized in columns named like `"<score>_original_target"`.

## How Pairing Is Determined

### Precursor-Level Pairs

- The spectral library carries a precomputed `:entrapment_pair_id` per precursor (e.g., peptide ion). That pairing is established during library preparation and ensures that each entrapment has a unique target counterpart.
- Use `add_entrap_pair_ids!(prec_results, library_precursors)` to attach the precomputed `:entrapment_pair_id` onto the results table via `:precursor_idx`.
- Once `:entrapment_pair_id` is present, `add_original_target_scores!` can compute `"<score>_original_target"` by locating the original target within the same pair (and file) and copying its score.

Notes
- Missing partners: If a paired original target does not exist in a file, the `"<score>_original_target"` value is set to `-1.0` for that row.
- Modifications: Utilities like `getModKey` normalize modification annotations for grouping during library construction; however, pairing at runtime simply consumes `:entrapment_pair_id` from the library.

### Protein-Level Pairs

- When `:entrapment_pair_id` is not present at the protein level, use `assign_protein_entrapment_pairs!` to create it from a table that has `:protein` and `:entrap_id`.
- Algorithm (per unique protein):
  - Split rows into original targets (`entrap_id == 0`) and entrapment groups (`entrap_id != 0`).
  - If a protein lacks either a target or an entrapment, it is skipped.
  - For each target row, assign a new pair ID and then round‑robin assign that same ID to exactly one row from each entrapment group. Round‑robin avoids bias when groups have unequal counts.
- Once assigned, propagate `:entrapment_pair_id` to results via `add_protein_entrap_pair_ids!` if needed (maps by `:protein`).
- Use `add_original_target_protein_scores!` to populate `"<score>_original_target"` at the protein level. The mapping keys on `(file, species, protein)` if `:species` is present, otherwise `(file, "", protein)`.
 - Use `add_original_target_protein_scores!` to populate `"<score>_original_target"` at the protein level. The mapping uses `(file, species, protein_key)` where `protein_key` is the first member of the protein group name (i.e., `split(protein,';')[1]`). If `:species` is absent, the key is `(file, "", protein_key)`.

Notes
- Duplicate targets: If a file/species/protein combination contains multiple targets, an error is thrown to prevent ambiguous mapping.
- Missing partners: If the original target is absent, the `"<score>_original_target"` value is set to `-1.0`.

## EFDR Methods

Two EFDR estimators are supported. Both operate on an ordering defined by `(-score, qval)` (descending score, then ascending q-value to break ties), and both enforce monotonicity of EFDR along that order.

### Combined EFDR

- Tracks cumulative counts:
  - `Nτ`: number of targets observed up to a cutoff.
  - `Nϵ`: number of entrapments observed up to a cutoff.
- EFDR estimate: `(Nϵ * (1 + 1/r)) / (Nϵ + Nτ)`, typically with `r = 1.0` for a balanced design.
- Implemented by constructing a `CombinedEFDR` method and calling `calculate_efdr`.

### Paired EFDR

- Uses original target scores to refine the numerator with additional terms that reflect within‑pair behavior:
  - `Nϵsτ`: entrapments that pass the entrapment cut but whose target is below the target cut (or missing).
  - `Nϵτs`: entrapments that exceed their target even when the target passes the cut.
- EFDR estimate: `(Nϵ + Nϵsτ + 2·Nϵτs) / (Nϵ + Nτ)`, clipped to [0, 1].
- Implemented by constructing a `PairedEFDR` method and calling `calculate_efdr`.
  - Exact computation can be O(n²); a sampling `stride` is supported for speed.
  - See also `calculate_efdr_fast(::PairedEFDR)` for an O(n log n + K) aggregated version.

## Global Best Rows Across Files

To create a single row per entity (precursor or protein), choose the best score across files:
- Precursors: `create_global_results_df(prec_results; score_col=:global_prob)` groups by `:precursor_idx`.
- Proteins: `create_global_protein_results_df(protein_results; score_col=:global_pg_score)` groups by `:species` (if present), `:protein`, and `:entrap_id` (if present), so that target and each entrapment group select their own best row.

Returned rows carry `:ms_file_idx = 0` or `:file_name = "global"` to indicate a cross-file selection.

## Practical Workflow

- Precursor EFDR
  1) Attach `:entrapment_pair_id` with `add_entrap_pair_ids!`.
  2) Compute original targets with `add_original_target_scores!` for desired score columns.
  3) Add EFDR columns via `add_efdr_columns!` (Combined always; Paired requires `"<score>_original_target"`).

- Protein EFDR
  1) Ensure `:entrapment_pair_id` exists: run `assign_protein_entrapment_pairs!` if necessary, or propagate from a library using `add_protein_entrap_pair_ids!`.
  2) Compute original targets with `add_original_target_protein_scores!`.
  3) Run `add_protein_efdr_columns!` to attach EFDR columns.

## Design Choices & Edge Cases

- Round‑robin pairing at the protein level ensures a fair distribution of pairs when entrapment groups have unequal number of rows.
- Missing partners are explicitly marked with `-1.0` in `"<score>_original_target"` columns; Paired EFDR treats these as target‑absent cases.
- EFDR curves are post‑processed to be monotonically nonincreasing along the `(-score, qval)` order.
- The ratio `r` defaults to 1.0 but can be tuned if the number of available target vs entrapment entities differs.

## References (Code Pointers)

- Pairing
  - `assign_protein_entrapment_pairs!` (protein): `src/core/protein_entrapment_pairing.jl`
  - `add_entrap_pair_ids!` (precursor): `src/core/entrapment_pairing.jl`
- Original targets
  - Precursors: `add_original_target_scores!` in `src/core/scoring.jl`
  - Proteins: `add_original_target_protein_scores!` in `src/core/protein_scoring.jl`
- EFDR methods
  - `CombinedEFDR`, `PairedEFDR`, and `add_efdr_columns!`: `src/core/efdr_methods.jl`
  - Fast Paired EFDR: `calculate_efdr_fast` in `src/core/paired_fast.jl`
