using Pkg
Pkg.add("CSV")

function Replace(file_path)
    p = CSV.read(file_path) |> Matrix
    println(p)

p = CSV.read("sample code/sample data/triblock-sequences.csv")
println(p)

# function KL_divergence(p, q)
#     kl = 0
#     for i in eachindex(p)
#         kl += p[i] * log2(p[i]/q[i])
#     end
#     return kl
# end    

# function JS_divergence(p,q)
#     p /= sum(p)
#     q /= sum(q)
#     m = (p + q) * 0.5
#     return 0.5 * KL_divergence(p, m) + 0.5 * KL_givergence(q,m)
# end

# Replace("sample code/sample data/triblock-sequences.csv")
