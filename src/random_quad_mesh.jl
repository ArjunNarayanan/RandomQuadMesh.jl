function random_quad_mesh(npts)
    pv = random_polygon(npts)
    return quad_mesh(pv)
end

function mwmatching_quad_mesh(boundary_points)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], Inf, "pQ")
    mesh = Mesh(p, t)
    q, t = match_tri2quad(mesh)
    quad_mesh = triquad_refine(mesh.p, q, t)
    return quad_mesh
end

function catmull_clark_quad_mesh(boundary_points)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], Inf, "pQ")
    q = zeros(Int, 4, 0)
    quad_mesh = triquad_refine(p, q, t)
    return quad_mesh
end

function quad_mesh(boundary_points; algorithm="matching")
    if algorithm == "matching"
        return mwmatching_quad_mesh(boundary_points)
    elseif algorithm == "catmull-clark"
        return catmull_clark_quad_mesh(boundary_points)
    else
        error("Expected algorithm = {matching, catmull-clark} but got algorithm = $algorithm")
    end
end

function tri_mesh(boundary_points)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], Inf, "pQ")
    mesh = Mesh(p, t)
end