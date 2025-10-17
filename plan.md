# Implementation Plan: Option 2 - Evidence-Based Protein Selection

## Overview

Implement evidence-based protein selection for entrapment pairing. Instead of using the first protein in a semicolon-delimited group, select the protein that appears most frequently in the supporting peptides' protein annotations.

## Key Changes

### 1. Add Helper Function: `get_best_protein_representative`
**Location:** `src/core/protein_scoring.jl`

**Purpose:** Given a protein group and its supporting peptides, return the single protein that appears most frequently in those peptides' protein_groups.

**Algorithm:**
1. Parse protein_group into candidate proteins
2. For each peptide in peptide_indices:
   - Look up peptide's :protein_groups in precursors_library
   - For each candidate protein, check if it appears in peptide's groups
   - Increment count for matching candidates
3. Return candidate with highest count

**Edge Cases:**
- Single-protein groups → Return immediately
- Empty peptide_indices → Fall back to first protein
- Missing peptide_indices → Fall back to first protein
- Invalid precursor_idx → Skip and continue
- Tie in counts → Use first occurrence order

### 2. Modify `add_original_target_protein_scores!`
**Location:** `src/core/protein_scoring.jl` (line ~34)

**Current signature:**
```julia
function add_original_target_protein_scores!(df::DataFrame,
                                              score_cols::Union{Symbol, Vector{Symbol}})
```

**New signature:**
```julia
function add_original_target_protein_scores!(df::DataFrame,
                                              score_cols::Union{Symbol, Vector{Symbol}};
                                              precursors_library::Union{Nothing, DataFrame}=nothing)
```

**Modified logic:**
```julia
# Determine protein key using evidence-based selection if available
protein_key = if precursors_library !== nothing && hasproperty(df, :peptides)
    # Use evidence-based selection
    peptide_indices = row.peptides
    get_best_protein_representative(row.protein, peptide_indices, precursors_library)
else
    # Fallback to first protein (current behavior)
    split(row.protein, ';')[1]
end
```

### 3. Update API Call Sites
**File:** `src/api.jl`

**Locations to update:**
- Line ~444: `add_original_target_protein_scores!` (per-file proteins, replicate analysis)
- Line ~460: `add_original_target_protein_scores!` (global proteins, replicate analysis)
- Line ~854: `add_original_target_protein_scores!` (per-file proteins, single analysis)
- Line ~874: `add_original_target_protein_scores!` (global proteins, single analysis)

**Change:**
```julia
# Before
add_original_target_protein_scores!(protein_results, [s for (s,_) in prot_perfile_pairs])

# After
add_original_target_protein_scores!(protein_results, [s for (s,_) in prot_perfile_pairs];
                                     precursors_library=library_precursors)
```

**Note:** For protein-only analyses (`run_protein_efdr_analysis`), library is not available. These will automatically fall back to first-protein method.

### 4. Add Comprehensive Tests
**Location:** `test/test_protein_entrapment.jl`

**Test coverage:**
- Basic selection by evidence count
- Single protein group
- Empty peptides (fallback)
- Missing peptides (fallback)
- Invalid indices (skip and continue)
- Tie in counts (first occurrence)
- Integration with `add_original_target_protein_scores!`
- Backward compatibility without library

## Implementation Checklist

- [ ] Read `src/core/protein_scoring.jl` to understand current implementation
- [ ] Add `get_best_protein_representative` function (~70 lines)
  - [ ] Handle single-protein groups
  - [ ] Handle missing/empty peptides
  - [ ] Initialize count dictionary
  - [ ] Query precursors_library for protein_groups
  - [ ] Count occurrences of each candidate
  - [ ] Return best protein with tie-breaking
- [ ] Modify `add_original_target_protein_scores!` signature
  - [ ] Add `precursors_library` keyword argument
  - [ ] Update protein_key selection logic
  - [ ] Maintain backward compatibility
- [ ] Update API call sites in `src/api.jl`
  - [ ] Line ~444 (per-file proteins, replicate)
  - [ ] Line ~460 (global proteins, replicate)
  - [ ] Line ~854 (per-file proteins, single)
  - [ ] Line ~874 (global proteins, single)
- [ ] Add comprehensive tests in `test/test_protein_entrapment.jl`
  - [ ] Test evidence-based selection
  - [ ] Test edge cases
  - [ ] Test integration
  - [ ] Test backward compatibility
- [ ] Run test suite to ensure no regressions
- [ ] Test with real data
  - [ ] Verify pairing improvements
  - [ ] Check performance overhead
  - [ ] Compare to current method

## Expected Benefits

✅ **Evidence-based:** Uses actual peptide-to-protein mapping
✅ **Robust pairing:** More likely to match targets with entrapments
✅ **Biologically meaningful:** Reflects peptide support
✅ **Backward compatible:** Falls back gracefully when library unavailable
✅ **Handles order variation:** Different protein orders can still pair correctly

## Performance Considerations

- Per protein group: O(num_peptides × avg_proteins_per_peptide)
- Typical case: 10 peptides × 2 proteins = 20 string comparisons
- Estimated overhead: ~0.1-1% of total analysis time

## Files to Modify

1. **src/core/protein_scoring.jl** (~100 lines total)
   - Add `get_best_protein_representative` function (~70 lines)
   - Modify `add_original_target_protein_scores!` (~10 lines changed)

2. **src/api.jl** (~8 lines changed)
   - Update 4 call sites to pass `precursors_library` parameter

3. **test/test_protein_entrapment.jl** (~140 lines added)
   - Add comprehensive test coverage

**Total estimated changes:** ~250 lines (code + tests)

## Open Questions

1. **Column name in library:** Does `precursors_table.arrow` use `:protein_groups`, `:protein`, or `:proteins`?
   - Current plan: Check for `:protein_groups` first, then `:protein`

2. **Peptides column format:** Is `:peptides` always `Vector{UInt32}`, or can it contain `missing` values?
   - Current plan: Handle both `Vector{UInt32}` and `Union{Missing, Vector{UInt32}}`

3. **Library availability:** For `run_protein_efdr_analysis`, should we add a library parameter?
   - Current plan: Use fallback method (keeps API simple)
