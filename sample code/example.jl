using DelimitedFiles
using CopolymerSequenceGeneration

println("Fetching feed ratio data from folder...")
feed_ratios = readdlm("sample code/sample data/triblock.csv", ',', Float64, skipstart=1)

println("Beginning sequence generation...")
results = run_seq(feed_ratios, 1000, 1, 1)

println("Exporting data to path.")
path = string("sample code/sample data/triblock-sequences.csv")
writedlm(path,  results.sequences, ',')
writedlm("sample code/sample data/triblock-length.csv", results.lengths, ',')