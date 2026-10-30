"""
Co-occurrence persistent homology over Z/2.

    include("cooccurrence_ph.jl")
    using .CooccurrencePH

    sentences = [["square", "odd"], ["square", "remainder"],
                 ["odd", "remainder"], ["odd", "square", "remainder"]]
    cooccurrence_persistence(sentences; maxdim=1)
    # [(0, 1, -1), (1, 3, 4)]

The input is a filtration given as a sequence of sets: step t adds the full simplex on the
t-th set (in the paper, the words of the t-th sentence). Vertex labels may be any sortable
values (strings, integers, ...). The output is the persistence diagram as sorted
`(dimension, birth, death)` triples for dimensions 0..maxdim. Steps are numbered from 1,
`death = -1` marks a class that never dies, and zero-length intervals are omitted.

Method (see Appendix A of the paper): dominated vertices and redundant steps are removed,
then a filtered acyclic matching (one apex vertex per step, chosen by earlier and total
frequency) reduces the filtration to its critical cells, whose boundary matrix is reduced
with clearing. Dimension 0 alone is computed with union-find. Standard library only.
"""
module CooccurrencePH

export cooccurrence_persistence

const Interval = Tuple{Int,Int,Int}

"""Relabel the vertices 1..nv (sorted by label); returns the sets as sorted Int32 vectors."""
function normalize_sets(sets)
    labels = sort!(unique!(collect(Iterators.flatten(sets))))
    index = Dict(v => Int32(i) for (i, v) in enumerate(labels))
    [sort!(unique!(Int32[index[v] for v in s])) for s in sets], length(labels)
end

@inline function find_root!(parent::Vector{Int32}, vertex::Int32)
    while parent[vertex] != vertex
        parent[vertex] = parent[parent[vertex]]
        vertex = parent[vertex]
    end
    vertex
end

"""Dimension 0 only: elder-rule union-find, one star per arriving set."""
function h0_homology(normalized::Vector{Vector{Int32}}, nv::Int)
    parent = zeros(Int32, nv)
    sizes = zeros(Int, nv)
    elder = Int32.(1:nv)
    births = zeros(Int, nv)
    out = Interval[]
    for (t, s) in enumerate(normalized)
        isempty(s) && continue
        for v in s
            if parent[v] == 0
                parent[v] = v
                sizes[v] = 1
                births[v] = t
            end
        end
        anchor = first(s)
        for v in Iterators.drop(s, 1)
            a, b = find_root!(parent, anchor), find_root!(parent, v)
            a == b && continue
            va, vb = elder[a], elder[b]
            older, younger = (births[va], va) <= (births[vb], vb) ? (va, vb) : (vb, va)
            births[younger] < t && push!(out, (0, births[younger], t))
            sizes[a] < sizes[b] && ((a, b) = (b, a))
            parent[b] = a
            sizes[a] += sizes[b]
            elder[a] = older
        end
    end
    for v in 1:nv
        parent[v] == v && push!(out, (0, births[elder[v]], -1))     # one class per component
    end
    sort!(out)
end

"""Remove dominated vertices and steps whose set is contained in an earlier set."""
function preprocess_sets(sets::Vector{Vector{Int32}}, nv::Int)
    occ = [Int[] for _ in 1:nv]
    for (t, s) in enumerate(sets), v in s
        push!(occ[v], t)
    end
    alive = trues(nv)
    for v in 1:nv
        candidates = copy(sets[first(occ[v])])
        for t in occ[v]
            filter!(w -> w in sets[t], candidates)
            length(candidates) <= 1 && break
        end
        # Strict degree / label orientation makes simultaneous deletions acyclic.
        for w in candidates
            if length(occ[w]) > length(occ[v]) || (length(occ[w]) == length(occ[v]) && w < v)
                alive[v] = false
                break
            end
        end
    end
    reduced = [Int32[v for v in s if alive[v]] for s in sets]
    prior = [Int[] for _ in 1:nv]
    for (t, s) in enumerate(reduced)
        isempty(s) && continue
        v = s[argmin([length(prior[w]) for w in s])]
        if any(i -> length(reduced[i]) >= length(s) && all(w -> w in reduced[i], s), prior[v])
            reduced[t] = Int32[]
        else
            for w in s
                push!(prior[w], t)
            end
        end
    end
    reduced
end

# Simplices are fixed-length tuples of vertex ids, padded with zeros.
@inline simplex_key(::Val{K}, vertices, n::Int) where K =
    ntuple(i -> i <= n ? Int32(vertices[i]) : Int32(0), Val(K))

@inline remove_vertex(key::NTuple{K,Int32}, at::Int) where K =
    ntuple(i -> i < at ? key[i] : i < K ? key[i+1] : Int32(0), Val(K))

@inline function add_vertex(key::NTuple{K,Int32}, n::Int, v::Int32) where K
    at = 1
    while at <= n && key[at] < v
        at += 1
    end
    ntuple(i -> i < at ? key[i] : i == at ? v : i <= n+1 ? key[i-1] : Int32(0), Val(K))
end

"""Call f on every n-subset of s, as a padded key."""
function each_combination(f, s::Vector{Int32}, n::Int, vk::Val{K}) where K
    m = length(s)
    m < n && return
    if n == 1
        for a in 1:m
            f(ntuple(i -> i == 1 ? s[a] : Int32(0), vk))
        end
    elseif n == 2
        for a in 1:m-1, b in a+1:m
            f(ntuple(i -> i == 1 ? s[a] : i == 2 ? s[b] : Int32(0), vk))
        end
    elseif n == 3
        for a in 1:m-2, b in a+1:m-1, c in b+1:m
            f(ntuple(i -> i == 1 ? s[a] : i == 2 ? s[b] : i == 3 ? s[c] : Int32(0), vk))
        end
    else
        buf = Vector{Int32}(undef, n)
        function grow(start, depth)
            if depth > n
                f(simplex_key(vk, buf, n))
                return
            end
            for j in start:m-(n-depth)
                buf[depth] = s[j]
                grow(j + 1, depth + 1)
            end
        end
        grow(1, 1)
    end
    nothing
end

"""Symmetric difference of two sorted vectors (addition of Z/2 chains), into out."""
function xor_sorted!(out::Vector{Int32}, a::Vector{Int32}, b::Vector{Int32})
    empty!(out)
    sizehint!(out, length(a) + length(b))
    i = j = 1
    @inbounds while i <= length(a) && j <= length(b)
        x, y = a[i], b[j]
        if x < y
            push!(out, x); i += 1
        elseif y < x
            push!(out, y); j += 1
        else
            i += 1; j += 1
        end
    end
    @inbounds while i <= length(a)
        push!(out, a[i]); i += 1
    end
    @inbounds while j <= length(b)
        push!(out, b[j]); j += 1
    end
    out
end

xor_sorted(a::Vector{Int32}, b::Vector{Int32}) =
    isempty(b) ? a : isempty(a) ? copy(b) : xor_sorted!(Int32[], a, b)

"""Birth of a simplex: the first step whose set contains all its vertices (bitset query)."""
function first_containing(key::NTuple{K,Int32}, incidence, tmax::Int) where K
    @inbounds for word in 1:cld(tmax, 64)
        bits = typemax(UInt64)
        for v in key
            v == 0 && break
            bits &= incidence[v][word]
            bits == 0 && break
        end
        bits != 0 && return 64 * (word - 1) + trailing_zeros(bits) + 1
    end
    typemax(Int)
end

struct MorseContext{K}
    birth::Dict{NTuple{K,Int32},Int}
    apex::Vector{Int32}
    incidence::Vector{Vector{UInt64}}
    nsets::Int
    index::Dict{NTuple{K,Int32},Int32}
    projection::Dict{NTuple{K,Int32},Vector{Int32}}
end

@inline birth!(ctx::MorseContext{K}, key::NTuple{K,Int32}) where K =
    get!(() -> first_containing(key, ctx.incidence, ctx.nsets), ctx.birth, key)

"""Projection of a simplex onto the critical cells (the recursion of Appendix A)."""
function project!(ctx::MorseContext{K}, key::NTuple{K,Int32}, n::Int) where K
    n == 0 && return Int32[]
    cached = get(ctx.projection, key, nothing)
    cached !== nothing && return cached
    t = birth!(ctx, key)
    p = ctx.apex[t]
    at = findfirst(==(p), key)
    result = Int32[]
    if at !== nothing
        if n == 1 || birth!(ctx, remove_vertex(key, at)) < t
            push!(result, ctx.index[key])
        end
    else
        for j in 1:n
            child = add_vertex(remove_vertex(key, j), n - 1, p)
            result = xor_sorted(result, project!(ctx, child, n))
        end
    end
    ctx.projection[key] = result
    result
end

"""Reduce the Morse boundary matrix of the critical cells (with clearing)."""
function reduce_critical_cells(cells::Vector{NTuple{K,Int32}}, times::Vector{Int},
                               dims::Vector{Int}, maxdim::Int, ctx) where K
    nc = length(cells)
    owner = zeros(Int32, nc)
    columns = [Int32[] for _ in 1:nc]
    killed = falses(nc)
    positive = falses(nc)
    out = Interval[]
    scratch = Int32[]
    for dim in maxdim+1:-1:0
        for j in eachindex(cells)
            (dims[j] != dim || killed[j]) && continue
            col = Int32[]
            if dim > 0
                for at in 1:dim+1
                    col = xor_sorted(col, project!(ctx, remove_vertex(cells[j], at), dim))
                end
            end
            while !isempty(col)
                low = last(col)
                prev = owner[low]
                if prev == 0
                    owner[low] = j
                    killed[low] = true
                    b, d = times[low], times[j]
                    dims[low] <= maxdim && b < d && push!(out, (dims[low], b, d))
                    break
                end
                xor_sorted!(scratch, col, columns[prev])
                col, scratch = scratch, col
            end
            if isempty(col)
                positive[j] = true
            else
                columns[j] = col
            end
        end
    end
    for j in eachindex(cells)
        positive[j] && !killed[j] && dims[j] <= maxdim && push!(out, (dims[j], times[j], -1))
    end
    sort!(out)
end

"""Critical cells of the filtered matching (apex p_t at each step), then reduction."""
function morse_homology(sets::Vector{Vector{Int32}}, nv::Int, maxdim::Int, ::Val{K}) where K
    vk = Val(K)
    ns = length(sets)
    apex = zeros(Int32, ns)
    frequency = zeros(Int, nv)            # total number of sets containing v
    past_frequency = zeros(Int, nv)       # number of earlier sets containing v
    for s in sets, v in s
        frequency[v] += 1
    end
    incidence = [zeros(UInt64, cld(ns, 64)) for _ in 1:nv]
    for (t, s) in enumerate(sets)
        isempty(s) && continue
        for v in s
            incidence[v][cld(t, 64)] |= UInt64(1) << ((t - 1) % 64)
        end
        # Apex rule: largest earlier count, then largest total count, then smaller label.
        apex[t] = s[argmax([(past_frequency[v], frequency[v]) for v in s])]
        for v in s
            past_frequency[v] += 1
        end
    end
    cells = NTuple{K,Int32}[]
    times = Int[]
    dims = Int[]
    for (t, s) in enumerate(sets)
        isempty(s) && continue
        p = apex[t]
        vertex = simplex_key(vk, (p,), 1)
        if first_containing(vertex, incidence, t) == t
            push!(cells, vertex); push!(times, t); push!(dims, 0)
        end
        base_vertices = Int32[v for v in s if v != p]
        for n in 1:min(K - 1, length(base_vertices))
            each_combination(base_vertices, n, vk) do base
                if first_containing(base, incidence, t) < t           # an old face ...
                    key = add_vertex(base, n, p)
                    if first_containing(key, incidence, t) == t       # ... plus the apex is new
                        push!(cells, key); push!(times, t); push!(dims, n)
                    end
                end
            end
        end
    end
    index = Dict(key => Int32(i) for (i, key) in enumerate(cells))
    ctx = MorseContext(Dict{NTuple{K,Int32},Int}(), apex, incidence, ns, index,
                       Dict{NTuple{K,Int32},Vector{Int32}}())
    reduce_critical_cells(cells, times, dims, maxdim, ctx)
end

"""
    cooccurrence_persistence(sets; maxdim=1) -> Vector{Tuple{Int,Int,Int}}

Persistence diagram, over Z/2 and in dimensions 0..maxdim, of the co-occurrence filtration
in which step t adds the full simplex on `sets[t]`. Returns sorted `(dimension, birth, death)`
triples; steps are numbered from 1, `death = -1` marks a class that never dies, and
zero-length intervals are omitted.
"""
function cooccurrence_persistence(sets; maxdim::Int=1)
    maxdim >= 0 || throw(ArgumentError("maxdim must be nonnegative"))
    all(isempty, sets) && return Interval[]
    normalized, nv = normalize_sets(sets)
    maxdim == 0 && return h0_homology(normalized, nv)
    reduced = preprocess_sets(normalized, nv)
    if maximum(length, reduced; init=0) <= 1
        # Only isolated vertices remain: one essential H0 class per vertex.
        seen = falses(nv)
        out = Interval[]
        for (t, s) in enumerate(reduced), v in s
            seen[v] || (push!(out, (0, t, -1)); seen[v] = true)
        end
        return sort!(out)
    end
    morse_homology(reduced, nv, maxdim, Val(maxdim + 2))
end

end # module
