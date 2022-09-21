function random_quad_mesh(npts)
    mesh = random_polygon_trimesh(npts)
    q, t = match_tri2quad(mesh)
    quad_mesh = triquad_refine(mesh.p, q, t)
    return quad_mesh
end