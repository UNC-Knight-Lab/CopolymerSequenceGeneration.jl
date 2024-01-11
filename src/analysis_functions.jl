function KL_divergence(p, q)
    kl = 0
    for i in eachindex(p)
        kl += p[i] * log2(p[i]/q[i])
    end
    return kl
end    

function JS_divergence(p,q)
    p /= sum(p)
    q /= sum(q)
    m = (p + q) * 0.5
    return 0.5 * KL_divergence(p, m) + 0.5 * KL_givergence(q,m)
end
