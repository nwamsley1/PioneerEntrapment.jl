using PioneerEntrapment
import PioneerEntrapment: calculate_efdr_fast

# Simple test case
e = Float64[0.95, 0.90, 0.85]
t = Float64[0.90, 0.92, 0.80]
labs = Int[0, 1, 0]  # 0=target, 1=entrapment
q = Float64[0.01, 0.02, 0.03]

println("Row | Entrap | Target | Label (0=targ, 1=entr)")
println("----+--------+--------+---------------------")
for i in 1:3
    label_str = labs[i] == 0 ? "target" : "entrap"
    println(lpad(i, 3), " | ", lpad(e[i], 6), " | ", lpad(t[i], 6), " | ", lpad(labs[i], 5), " ($label_str)")
end
println()

method = PairedEFDR(e, t, labs, q, 1.0)

# Manually compute what fast method should give
cuts = sort(unique(vcat(t, e)))
println("Cut points: ", cuts)
println()

for (j, cut) in enumerate(cuts)
    # Count at this threshold - ONLY count by label!
    T = sum((labs .== 0) .& (e .>= cut))  # TARGET rows with score >= cut
    E = sum((labs .== 1) .& (e .>= cut))  # ENTRAPMENT rows with score >= cut

    EST = 0  # entrap rows: score >= cut AND target_score < cut
    ETS = 0  # entrap rows: score >= cut AND target_score >= cut AND score > target_score

    for i in 1:3
        if labs[i] == 1  # Only check entrapment rows
            if e[i] >= cut
                if t[i] < cut
                    EST += 1
                elseif t[i] >= cut && e[i] > t[i]
                    ETS += 1
                end
            end
        end
    end

    efdr = (E + EST + 2*ETS) / (T + E)
    println("Cut $j: s=$cut")
    println("  T=$T (target rows >= cut), E=$E (entrap rows >= cut)")
    println("  EST=$EST, ETS=$ETS")
    println("  EFDR = ($E + $EST + 2*$ETS) / ($T + $E) = $efdr")
    println()
end

# Now see how fast method maps to rows
println("=" ^ 60)
println("Fast method row mapping:")
println("=" ^ 60)
fast_all = calculate_efdr_fast(method; cuts_mode=:all)
for i in 1:3
    ie = searchsortedlast(cuts, e[i])
    println("Row $i: e=$(e[i]) → cut index $ie (cut=$(ie > 0 ? cuts[ie] : "N/A")) → EFDR=$(fast_all[i])")
end
println()

standard = calculate_efdr(method; stride=1)
println("Standard EFDR: ", standard)
println("Fast EFDR:     ", fast_all)
println("Difference:    ", abs.(standard .- fast_all))
