using Revise
using RandomQuadMesh
using PlotQuadMesh

RM = RandomQuadMesh
PQ = PlotQuadMesh


# using Random
# Random.seed!(123)


boundary_points = RM.random_polygon(10)
quad_mesh = RM.quad_mesh(boundary_points, algorithm = "catmull-clark")

PQ.plot_mesh(quad_mesh.p, quad_mesh.t)