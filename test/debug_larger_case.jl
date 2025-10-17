using PioneerEntrapment
using Random
import PioneerEntrapment: calculate_efdr_fast

# Test from the suite that's failing
Random.seed!(42)
e = Float64[0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50]
t = Float64[0.90, 0.92, 0.80, 0.82, 0.70, 0.72, 0.50, 0.58, 0.45, 0.40]
labs = Int[0, 1, 0, 1, 0, 1, 1, 0, 1, 0]
q = Float64[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]

println("Row | Entrap | Target | Label | Q-val")
println("----+--------+--------+-------+------")
for i in 1:length(e)
    label_str = labs[i] == 0 ? "T" : "E"
    println(lpad(i, 3), " | ", lpad(e[i], 6), " | ", lpad(t[i], 6), " | ", lpad(label_str, 5), " | ", lpad(q[i], 5))
end
println()

method = PairedEFDR(e, t, labs, q, 1.0)
standard = calculate_efdr(method; stride=1)
fast_all = calculate_efdr_fast(method; cuts_mode=:all)

println("Standard EFDR: ", standard)
println("Fast EFDR:     ", fast_all)
println("Difference:    ", abs.(standard .- fast_all))
println("Max diff:      ", maximum(abs.(standard .- fast_all)))
println()

# Show where they differ
for i in 1:length(e)
    diff = abs(standard[i] - fast_all[i])
    if diff > 1e-10
        println("Row $i: standard=$(standard[i]), fast=$(fast_all[i]), diff=$diff")
    end
end
