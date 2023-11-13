using Revise
using RandomQuadMesh
using PlotQuadMesh
RM = RandomQuadMesh
PQ = PlotQuadMesh

using Random
Random.seed!(123)


allow_vertex_insert = true
hmax = Inf
algorithm = "catmull-clark"
# algorithm = "matching"

poly = RM.random_polygon(10)
mesh = RM.quad_mesh(
    poly,
    algorithm=algorithm,
    hmax=hmax,
    allow_vertex_insert=allow_vertex_insert
)


fig, ax = PQ.plot_mesh(mesh.p, mesh.t)
fig

# mesh = RM.random_polygon_trimesh(10)
# q, t = RM.match_tri2quad(mesh)
# quad_mesh = RM.triquad_refine(mesh.p, q, t)

# using QuadMeshGame
# QM = QuadMeshGame
# qmesh = QM.QuadMesh(quad_mesh.p, quad_mesh.t, quad_mesh.t2t, quad_mesh.t2n)