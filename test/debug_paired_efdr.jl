using PioneerEntrapment
import PioneerEntrapment: calculate_efdr_fast

# Simple test case
e = Float64[0.95, 0.90, 0.85]
t = Float64[0.90, 0.92, 0.80]
labs = Int[0, 1, 0]
q = Float64[0.01, 0.02, 0.03]

println("Entrap scores (e): ", e)
println("Target scores (t): ", t)
println("Labels: ", labs)
println("Q-values: ", q)
println()

method = PairedEFDR(e, t, labs, q, 1.0)

println("=" ^ 60)
println("Testing with cuts_mode = :targets_only (default)")
println("=" ^ 60)
standard = calculate_efdr(method; stride=1)
fast_targets = calculate_efdr_fast(method; cuts_mode=:targets_only)

println("Standard EFDR:       ", standard)
println("Fast EFDR (targets): ", fast_targets)
println("Difference:          ", abs.(standard .- fast_targets))
println("Max diff:            ", maximum(abs.(standard .- fast_targets)))
println()

println("=" ^ 60)
println("Testing with cuts_mode = :all")
println("=" ^ 60)
fast_all = calculate_efdr_fast(method; cuts_mode=:all)

println("Standard EFDR:    ", standard)
println("Fast EFDR (all):  ", fast_all)
println("Difference:       ", abs.(standard .- fast_all))
println("Max diff:         ", maximum(abs.(standard .- fast_all)))
