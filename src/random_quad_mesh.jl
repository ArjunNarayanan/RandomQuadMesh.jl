function random_quad_mesh(npts)
    pv = random_polygon(npts)
    return quad_mesh(pv)
end

function quad_mesh(boundary_points)
    p, t = polytrimesh([boundary_points[:,[1:end;1]]], [], Inf, "pQ")
    mesh = Mesh(p, t)
    q, t = match_tri2quad(mesh)
    quad_mesh = triquad_refine(mesh.p, q, t)
    return quad_mesh
end