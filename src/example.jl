using DelimitedFiles

mutable struct System_features
    sequences::Array{Int64}
    lengths::Array{Int64}
    max_DP::Int64
    n_chains::Int64
    r_A::Float64
    r_B::Float64
end

feed_ratios = readdlm("pkg_files/sample data/triblock.csv", ',', Float64, skipstart=1)
num_blocks = length(feed_ratios[:,1])

System_features.max_DP = round(sum(feed_ratios)) + 50
System_features.n_chains = 1000
System_features.sequences = zeros(System_features.n_chains, System_features.max_DP)
System_features.lengths = ones(Int64, System_features.n_chains)
System_features.r_A = 1
System_features.r_B = 1

p = Progress(num_blocks; 1)

for block in 1:num_blocks
    A_mmol = feed_ratios[block, 1]
    B_mmol = feed_ratios[block, 2]

    
end

# path = string("/Users/suprajachittari/Documents/peter/sequence/121423/statistical-full.csv")
# writedlm(path,  all_seqs, ',')

# writedlm( "/Users/suprajachittari/Documents/peter/sequence/121423/statistical-lengths.csv",  position, ',')