using Revise
using Test
using RandomQuadMesh
RM = RandomQuadMesh

using Random
Random.seed!(123)


mesh = RM.random_polygon_trimesh(10)
q, t = RM.match_tri2quad(mesh)

quad_mesh = RM.triquad_refine(mesh.p, q, t)


using 