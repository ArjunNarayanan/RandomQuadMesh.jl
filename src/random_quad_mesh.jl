function random_quad_mesh(npts)
    mesh = RM.random_polygon_trimesh(npts)
    q, t = RM.match_tri2quad(mesh)
    quad_mesh = RM.triquad_refine(mesh.p, q, t)
    return quad_mesh
end