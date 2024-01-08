const volume = 1

"""sequence of all chains as well as length of each chain"""
mutable struct System_features
    max_DP::Int64
    n_chains::Int64
    sequences::Array{Int64}
    lengths::Array{Int64}
    r_A::Float64
    r_B::Float64
end

function first_monomer(A::Float64, B::Float64)
    m_A = A / volume
    m_B = B / volume

    if m_A > m_B # choose monomer A
        return 1 #0
    else # choose monomer B
        return 2 #1
    end
end

function growth(A::Float64, B::Float64, sequence::Array{Int64}, live_position::Int64, r_A::Float64, r_B::Float64)
    last_monomer = sequence[live_position - 1]

    rates_add = zeros(2,3) # first column is normalized rate, second is addition, third is sum

    if last_monomer == 0 # last monomer is A 
        k_same = 1 # rate constant for adding more A
        k_new = r_A * k_same # rate constant for adding B

        norm_rates = (k_same * (A / volume)) + (k_new * (B / volume))
        rates_add[1,1] = (k_same * (A / volume)) / norm_rates
        rates_add[1,2] = 1
        rates_add[1,3] = sum(rates_add[:,1])
        rates_add[2,1] = (k_new * (B / volume)) / norm_rates
        rates_add[2,2] = 2
        rates_add[2,3] = sum(rates_add[:,1])

    else 
        k_same = 1 # rate constant for adding more B
        k_new = r_B * k_same # rate constant for adding A

        norm_rates = (k_same * (B / volume)) + (k_new * (A / volume))
        rates_add[1,1] = (k_same * (B / volume)) / norm_rates
        rates_add[1,2] = 2
        rates_add[1,3] = sum(rates_add[:,1])
        rates_add[2,1] = (k_new * (A / volume)) / norm_rates
        rates_add[2,2] = 1
        rates_add[2,3] = sum(rates_add[:,1])
    end
    # println(rates_add)
    u = rand()*rates_add[2,3]

    if rates_add[1,3] < u <= rates_add[2,3]
        return Int(rates_add[2,2])
    elseif u <= rates_add[1,3]
        return Int(rates_add[1,2])
    else
        println("Error performing move, random number u is ", u)
    end
end

function update_system(A::Float64, B::Float64, move::Int64, live_position::Int64, sequence::Array{Int64}, CTA_mmol::Float64)
    change = CTA_mmol / volume

    if move == 1 #0 # move is adding A 
        if A - change < 0
            return sequence, A, B, live_position
        else
            sequence[live_position] = 1 #0
            A -= change
            live_position += 1
            return sequence, A, B, live_position
        end
    else
        if B - change < 0
            return sequence, A, B, live_position
        else
            sequence[live_position] = 2 #1
            B -= change
            live_position += 1
            return sequence, A, B, live_position
        end
    end
end


function run_block(System_features, A_mmol::Float64, B_mmol::Float64)
    #all_seqs::Array{Int64}, position::Array{Int64}
    all_seqs = System_features.sequences
    position = System_features.lengths
    max_DP = System_features.max_DP
    n_chains = System_features.n_chains
    r_A = System_features.r_A
    r_B = System_features.r_B

    CTA_mmol = 1 / n_chains

    attempt = 1

    while minimum(position) <= max_DP

        if attempt > n_chains * max_DP + 10000
            break
        end

        chain = rand(1:n_chains)

        if position[chain] == 1 # special case for first monomer
            new = first_monomer(A_mmol, B_mmol)
            all_seqs[chain, :], A_mmol, B_mmol, position[chain] = update_system(A_mmol, B_mmol, new, position[chain], all_seqs[chain, :], CTA_mmol)
        else
            new = growth(A_mmol, B_mmol, all_seqs[chain, :], position[chain], r_A, r_B)
            all_seqs[chain, :], A_mmol, B_mmol, position[chain] = update_system(A_mmol, B_mmol, new, position[chain], all_seqs[chain, :], CTA_mmol)
        end

        # println(A_mmol, B_mmol)

        if A_mmol < 0 || B_mmol < 0
            println("ERROR: negative concentrations obtained.")
            break
        end

        attempt += 1
    end

    println("Simulation complete. Final concentrations are ", A_mmol, " and ", B_mmol)

    System_features.sequences = all_seqs
    System_features.lengths = position

    return System_features
end

function run_seq(feed_ratios, n_chains, r_A, r_B)
    num_blocks = length(feed_ratios[:,1])
    max_DP = round(Int64, sum(feed_ratios)) + 50
    sequences = zeros(n_chains, max_DP)
    lengths = ones(Int64, n_chains)
    system = System_features(max_DP, n_chains, sequences, lengths, r_A, r_B)

    p = Progress(num_blocks; dt=1.0)

    for block in 1:num_blocks
        A_mmol = feed_ratios[block, 1]
        B_mmol = feed_ratios[block, 2]

        system = run_block(system, A_mmol, B_mmol)
        
        next!(p)
    end

    return system
end