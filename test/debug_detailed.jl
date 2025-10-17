using PioneerEntrapment
import PioneerEntrapment: calculate_efdr_fast

# Simple test case
e = Float64[0.95, 0.90, 0.85]
t = Float64[0.90, 0.92, 0.80]
labs = Int[0, 1, 0]
q = Float64[0.01, 0.02, 0.03]

println("Row | Entrap | Target | Label | Q-val")
println("----+--------+--------+-------+------")
for i in 1:length(e)
    println(lpad(i, 3), " | ", lpad(e[i], 6), " | ", lpad(t[i], 6), " | ", lpad(labs[i], 5), " | ", lpad(q[i], 5))
end
println()

method = PairedEFDR(e, t, labs, q, 1.0)

# Manually trace the standard method logic
println("=" ^ 60)
println("Standard Method (stride=1) - Tracing Logic")
println("=" ^ 60)
sort_indices = sortperm(collect(zip(-method.score, method.qval)))
println("Sorted indices (by -score, qval): ", sort_indices)
println()

for k in sort_indices
    s = method.score[k]
    println("Position k=$k, threshold s=$s (from row $k)")

    # Count logic from standard method
    Nτ = 0  # targets >= s
    Nϵ = 0  # entrapments >= s
    Nϵsτ = 0  # entrap >= s with target < s
    Nϵτs = 0  # entrap > target with target >= s

    for j in sort_indices[1:findfirst(==(k), sort_indices)]
        is_target = method.entrapment_label[j] == 0
        e_j = method.score[j]
        t_j = method.original_target_score[j]

        if is_target
            Nτ += 1
        else
            Nϵ += 1
            if (e_j >= s) && (t_j < s)
                Nϵsτ += 1
            elseif (e_j > t_j) && (t_j >= s)
                Nϵτs += 1
            end
        end
    end

    efdr = (Nϵ + Nϵsτ + 2*Nϵτs) / (Nϵ + Nτ)
    println("  Nτ=$Nτ, Nϵ=$Nϵ, Nϵsτ=$Nϵsτ, Nϵτs=$Nϵτs")
    println("  EFDR = ($Nϵ + $Nϵsτ + 2*$Nϵτs) / ($Nϵ + $Nτ) = $efdr")
    println()
end

standard = calculate_efdr(method; stride=1)
println("Final standard EFDR: ", standard)
println()

println("=" ^ 60)
println("Fast Method - What cut points are used?")
println("=" ^ 60)
cuts_all = sort(unique(vcat(t, e)))
println("Cuts (all mode): ", cuts_all)
println()
