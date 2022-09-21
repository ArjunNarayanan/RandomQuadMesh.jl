function py_maxWeightMatching(g)
    current = pwd()
    py"""
    import sys
    import os
    path = os.path.join(os.getcwd(), "src")
    sys.path.append(path)
    from mwmatching import *
    """
    out = py"maxWeightMatching"(g) .+ 1
end

function dual_mesh_edges(m::Mesh)
    g = []
    for ie = 1:size(m.t2t,1), it = 1:size(m.t2t,2)
        if m.t2t[ie,it] > it
            push!(g, ( it-1, m.t2t[ie,it]-1, 1 ))
        end
    end
    g
end

function match_tri2quad(m::Mesh)
    # Maximum matching on dual mesh
    g = dual_mesh_edges(m)
    out = py_maxWeightMatching(g)

    is_triangular(m) || throw("Input must be triangular mesh")
    
    # Form quads from matched triangles
    map = edgemap(m)
    q = []
    for it = 1:length(out)
        jt = out[it]
        if jt > it
            j = findfirst(m.t2t[:,it] .== jt)
            k = m.t2n[j,it]
            prev,next = [3,1,2],[2,3,1]
            cq = [ m.t[map[:,next[j]],it];
                   m.t[map[:,next[k]],jt] ]
            push!(q,cq)
        end
    end
    q = hcat(q...)

    # Unmatched triangles
    t = m.t[:,out.==0]

    q,t
end


function triquad_refine(p,q,t)
    # Refine hybrid tri/quad mesh into all-quad mesh
    
    qmap,tmap = edgemap(Quadrilateral()), edgemap(Triangle())
    nq,nt = size(q,2), size(t,2)

    dd = Dict{Tuple{Int,Int}, Int}()
    sizehint!(dd, nt*3 + nt*4)

    function el2edge(it, ie, el, map)
        e1 = el[map[1,ie],it]
        e2 = el[map[2,ie],it]
        e = ( min(e1,e2), max(e1,e2) )
    end

    # Edge mid-points (no duplication)
    emid = Vector{Float64}[]
    function addemid(e)
        if !haskey(dd,e)
            newmid = (p[:,e[1]] + p[:,e[2]]) / 2
            push!(emid, newmid)
            dd[e] = length(emid)
        end
    end
    for it = 1:nt, ie = 1:size(tmap,2)
        addemid(el2edge(it, ie, t, tmap))
    end
    for iq = 1:nq, ie = 1:size(qmap,2)
        addemid(el2edge(iq, ie, q, qmap))
    end

    # Element mid-points
    tmid = [ sum(p[:,ct],dims=2)/3 for ct in eachcol(t) ]
    qmid = [ sum(p[:,cq],dims=2)/4 for cq in eachcol(q) ]

    # Form new quads
    qadd = Vector{Int}[]
    np = size(p,2)
    ne = length(emid)

    # ... from triangles
    for it = 1:nt
        nds = [ dd[el2edge(it, ie, t, tmap)] + np for ie = 1:3 ]
        mid = np + ne + it
        newqs = [ [ t[1,it], nds[3], mid, nds[2] ],
                  [ t[2,it], nds[1], mid, nds[3] ],
                  [ t[3,it], nds[2], mid, nds[1] ] ]
        append!(qadd, newqs)
    end

    # ... from quads
    for iq = 1:nq
        nds = [ dd[el2edge(iq, ie, q, qmap)] + np for ie = 1:4 ]
        mid = np + ne + nt + iq
        newqs = [ [ q[1,iq], nds[1], mid, nds[4] ],
                  [ q[2,iq], nds[2], mid, nds[1] ],
                  [ q[3,iq], nds[3], mid, nds[2] ],
                  [ q[4,iq], nds[4], mid, nds[3] ] ]
        append!(qadd, newqs)
    end

    # Final mesh
    newp = hcat(p, emid..., tmid..., qmid...)
    newq = hcat(qadd...)
    Mesh(newp, newq)
end
