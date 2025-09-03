# Fast Paired EFDR: How it Works and Behavioral Equivalence

This note explains an O(n log n + K) algorithm for the paired entrapment EFDR that matches the current O(n^2) implementation. It summarizes how counts are accumulated across score cutoffs and why the resulting EFDR curve is identical (after the same monotone post-processing) to the quadratic method.

## Background

For each entity, we have:
- `e`: entrap score
- `t`: original target score (may be missing/NA)
- `label`: 0 for original target rows, 1 for entrap rows (we only need `t` and `e` for counting)

At a threshold `s`, the paired EFDR numerator is:
- `n_e(s)`: entraps with `e ≥ s`
- `n_est(s)`: entraps with `e ≥ s` and `(t < s or t is NA)`
- `n_ets(s)`: both discovered and `e > t` among pairs with `e ≥ s` and `t ≥ s`

The denominator is `n_targets(s) + n_e(s)` where `n_targets(s) = #(t ≥ s)`.

EFDR(s) = (n_e + n_est + 2*n_ets) / (n_targets + n_e)

The O(n^2) reference scans the sorted list and recomputes these counts at every step. The fast method pre-aggregates over a discrete set of ascending cutoffs and uses difference arrays to recover counts for all `s` in a single pass.

## Cutoffs and Indexing

Choose ascending cutoffs `cuts[1..K]`. Two sensible choices:
- targets_only: `cuts = sort(unique(t[!missing]))`
- all: `cuts = sort(unique(t ∪ e))`

For any value `x`, define `idx(x) = max{ j | cuts[j] ≤ x }` (0 if `x < cuts[1]`). With this, `x ≥ cuts[j]` iff `idx(x) ≥ j`.

We compute index vectors:
- `idx_t_all[i] = idx(t[i])` or missing
- `idx_e_all[i] = idx(e[i])` or missing

## Suffix Sums for n_targets and n_e

Let `tab(idx_vec)` produce a length-K histogram of 1..K indices (ignore 0/missing). Then:
- `T_asc = suffix_sum(tab(idx_t))` equals `n_targets(s)` for each ascending cutoff.
- `E_asc = suffix_sum(tab(idx_e))` equals `n_e(s)` for each ascending cutoff.

Suffix sum converts counts-at-exact-cut into counts-above-or-equal-to-cut.

## n_ets(s): entrap outranks target

Among rows where both scores are present and `e > t`, let `m_i = min(idx_t[i], idx_e[i])`. Each such pair contributes +1 to `n_ets(s)` for all cuts `≤ m_i`. Thus:
- histogram the `m_i` (ignoring 0)
- suffix sum that histogram to obtain `n_ets(s)` at each ascending cutoff.

## n_est(s): entrap ≥ s and (target < s or target is NA)

Think in terms of ranges on the ascending index `j` (cut position):
- If `t` is NA and `idx_e = k > 0`, the row contributes +1 for cuts `j = 1..k`.
- If both present with `idx_t = a`, `idx_e = b`, and `a+1 ≤ j ≤ b`, it contributes +1 for those `j`.

We accumulate these via a difference array `diff[1..K+1]`:
- For each segment `[start..end]`, do `diff[start] += 1; diff[end+1] -= 1`.
- Then `cumsum(diff)[1:K]` yields `EST_asc`.

This turns many range updates into O(K) work overall.

## EFDR and Mapping Back

With `T_asc, E_asc, EST_asc, ETS_asc` for ascending cuts, we compute `efdr_asc`. To report EFDR per row (for plotting and API compatibility), we map each item to its entrap index `idx_e_all[i]` and take `efdr_asc[idx_e_all[i]]`. Finally, we apply the same monotone correction used by the quadratic code over the sort order (by `(-score, qval)`).

## Should the Fast Method Match the Quadratic Method?

Yes, provided that:
- The same notion of threshold `s` is used (we evaluate at cutoffs consistent with the quadratic’s traversal), and
- The same monotone post-processing (reverse cumulative minimum over the traversal order) is applied.

Intuitively, both methods evaluate the same counts but the fast method reorganizes the computation: rather than recomputing counts at each step, it builds them once over all cutoffs via histograms, suffix sums, and a difference array. The EFDR curve (as a function of `s`) is identical. Ties and traversal order are reconciled by applying the same monotonicity correction.

## Edge Cases

- Missing target score: handled in `n_est` by giving full range `[1..idx_e]`.
- `idx = 0` (score below the first cutoff): contributes nothing to counts at positive cuts.
- Empty cut list: return all-zero EFDR.

## Complexity

- Building `cuts`: `O(n log n)` for sorting unique scores.
- Tabulation/suffix sums/difference array: `O(n + K)`
- Mapping per-row EFDR and monotone pass: `O(n)`

This replaces the O(n^2) reference with `O(n log n + K)`, where `K` is the number of unique cutoffs.

