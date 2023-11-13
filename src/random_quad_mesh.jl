function random_quad_mesh(npts)
    pv = random_polygon(npts)
    return quad_mesh(pv)
end

function mwmatching_quad_mesh(
    boundary_points,
    hmax,
    triangle_alg
)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], hmax, triangle_alg)
    mesh = Mesh(p, t)
    q, t = match_tri2quad(mesh)
    quad_mesh = triquad_refine(mesh.p, q, t)
    return quad_mesh
end

function catmull_clark_quad_mesh(
    boundary_points,
    hmax,
    triangle_alg
)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], hmax, triangle_alg)
    q = zeros(Int, 4, 0)
    quad_mesh = triquad_refine(p, q, t)
    return quad_mesh
end

function quad_mesh(boundary_points; algorithm="matching", hmax = Inf, allow_vertex_insert = false)
    triangle_alg = get_triangle_command(allow_vertex_insert)

    if algorithm == "matching"
        return mwmatching_quad_mesh(boundary_points, hmax, triangle_alg)
    elseif algorithm == "catmull-clark"
        return catmull_clark_quad_mesh(boundary_points, hmax, triangle_alg)
    else
        error("Expected algorithm = {matching, catmull-clark} but got algorithm = $algorithm")
    end
end

function get_triangle_command(allow_vertex_insert)
    if allow_vertex_insert
        return "puq28.6Q"
    else
        return "pQ"
    end
end

function tri_mesh(boundary_points; hmax = Inf, allow_vertex_insert = true)
    cmd = get_triangle_command(allow_vertex_insert)
    
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], hmax, cmd)
    return Mesh(p, t)
end