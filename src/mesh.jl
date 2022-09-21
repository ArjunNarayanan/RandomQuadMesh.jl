abstract type ElementType end
struct Triangle <: ElementType end
struct Quadrilateral <: ElementType end

name(::Triangle) = "triangular"
name(::Quadrilateral) = "quadrilateral"

struct Mesh
    p::Matrix{Float64}
    t::Matrix{Int64}
    eltype::ElementType
    t2t::Matrix{Int64}
    t2n::Matrix{Int64}
end

mkmatrix(x) = reshape(x, size(x,1), size(x,2)) # For single-column matrices
function Mesh(p,t,eltype)
    t2t,t2n = mkt2t(t,eltype)
    Mesh(mkmatrix(p), mkmatrix(t), eltype, mkmatrix(t2t), mkmatrix(t2n))
end

function Mesh(p, t)
    if size(t,1) == 3
        eltype = Triangle()
    elseif size(t,1) == 4
        eltype = Quadrilateral()
    else
        error("Unknown element type")
    end
    Mesh(p, t, eltype)
end

edgemap(::Triangle) =
    [2 3 1;
     3 1 2]
edgemap(::Quadrilateral) =
    [1 2 3 4;
     2 3 4 1]
edgemap(m::Mesh) = edgemap(m.eltype)

is_triangular(m::Mesh) = isa(m.eltype, Triangle)
is_quadrilateral(m::Mesh) = isa(m.eltype, Quadrilateral)

function Base.show(io::IO, m::Mesh)
    println(io, "Mesh:")
    println(io, "  $(size(m.p,2)) nodes")
    println(io, "  $(size(m.t,2)) $(name(m.eltype)) elements")
end

function mkt2t(t, eltype)
    map = edgemap(eltype)
    ne,nt = size(map,2), size(t,2)

    t2t = zeros(Int, ne, nt)
    t2n = zeros(Int, ne, nt)
    dd = Dict{Tuple{Int,Int}, Tuple{Int,Int}}()
    sizehint!(dd, nt*ne)
    for it = 1:nt
        for ie = 1:ne
            e1 = t[map[1,ie],it]
            e2 = t[map[2,ie],it]
            e = ( min(e1,e2), max(e1,e2) )
            if haskey(dd,e)
                nb = pop!(dd,e)
                t2t[ie,it] = nb[1]
                t2n[ie,it] = nb[2]
                t2t[nb[2],nb[1]] = it
                t2n[nb[2],nb[1]] = ie
            else
                dd[e] = (it,ie)
            end
        end
    end
    t2t,t2n
end

function boundedges(m::Mesh)
    map = edgemap(m)
    tedges = ( sort(ct[edge]) for edge = eachcol(map), ct = eachcol(m.t) )

    # Insert edges into set - if an edge occurs more than once, not on boundary
    ee = Set{Vector{Int}}()
    for e in tedges
        if e in ee
            delete!(ee, e)
        else
            push!(ee, e)
        end
    end
    collect(ee)
end

function boundary_interior_nodes(m::Mesh)
    bedges = boundedges(m)
    bnodes = unique(vcat(bedges...))
    inodes = setdiff(1:size(m.p,2), bnodes)
    bnodes,inodes
end

function uniref(m::Mesh, nref)
    m = deepcopy(m)
    for i = 1:nref
        m = uniref(m)
    end
    m
end

function uniref(m::Mesh)
    tedges = sort.([ ct[edge] for edge = eachcol(edgemap(m)), ct = eachcol(m.t) ])
    edges = unique(tedges)
    map = indexin(tedges, edges)

    np = size(m.p,2)
    nt = size(m.t,2)
    nedges = length(edges)

    pmid = [ sum(m.p[dim,edge])/2 for dim = 1:2, edge in edges ]
    mapmid = map .+ np
    newp = hcat(m.p, pmid)
    
    if is_triangular(m)
        prev = [3,1,2]
        next = [2,3,1]
        ts = [ [m.t[i,it],mapmid[prev[i],it],mapmid[next[i],it]] for i = 1:3, it = 1:nt ]
        ts2 = [ mapmid[:,it] for it = 1:nt ]
        ts = vcat(ts[:], ts2[:])
    elseif is_quadrilateral(m)
        pc = [ sum(m.p[dim,el])/4 for dim = 1:2, el = eachcol(m.t) ]
        mapc = np + nedges .+ (1:nt)
        newp = hcat(newp, pc)
        prev = [4,1,2,3]
        ts = [ [m.t[i,:] mapmid[i,:] mapc mapmid[prev[i],:]]' for i = 1:4 ]
    end
    newt = hcat(ts...)

    return Mesh(newp,newt)
end

function quad_area(quad)
    x1,y1,x2,y2,x3,y3,x4,y4 = quad

    AC = ( x3-x1, y3-y1 )
    BD = ( x4-x2, y4-y2 )

    A = 0.5*abs(AC[1]*BD[2] - AC[2]*BD[1])
end

function tri_qual(tri)
    x1,y1,x2,y2,x3,y3 = tri

    a = sqrt((x2-x1)^2 + (y2-y1)^2)
    b = sqrt((x3-x1)^2 + (y3-y1)^2)
    c = sqrt((x3-x2)^2 + (y3-y2)^2)

    r = 0.5*sqrt((b+c-a) * (c+a-b) * (a+b-c) / (a+b+c))
    R = a*b*c / sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))
    q = 2*r/R
end

function quad_kappa(quad)
    x1,y1,x2,y2,x3,y3,x4,y4 = quad

    a2 = (x2-x1)^2 + (y2-y1)^2
    b2 = (x3-x2)^2 + (y3-y2)^2
    c2 = (x4-x3)^2 + (y4-y3)^2
    d2 = (x1-x4)^2 + (y1-y4)^2
    a,b,c,d = sqrt(a2),sqrt(b2),sqrt(c2),sqrt(d2)

    sal = ((x1-x4)*(y2-y1) - (y1-y4)*(x2-x1)) / (d * a)
    sbe = ((x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)) / (a * b)
    sga = ((x3-x2)*(y4-y3) - (y3-y2)*(x4-x3)) / (b * c)
    sde = ((x4-x3)*(y1-y4) - (y4-y3)*(x1-x4)) / (c * d)

    ka = (d2 + a2) / (d * a * sal)
    kb = (a2 + b2) / (a * b * sbe)
    kc = (b2 + c2) / (b * c * sga)
    kd = (c2 + d2) / (c * d * sde)

    if ka<0 || kb<0 || kc<0 || kd<0
        return 0.0
    else
        return 4.0 / sqrt(ka^2 + kb^2 + kc^2 + kd^2)
    end
end

function elemqual(m::Mesh, fqual=nothing)
    if fqual == nothing
        if is_triangular(m)
            fqual = tri_qual
        elseif is_quadrilateral(m)
            fqual = quad_kappa
        end
    end

    [ fqual(m.p[:,el]) for el in eachcol(m.t) ]
end
