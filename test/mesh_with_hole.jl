using RandomQuadMesh
using MeshPlotter
using PlotQuadMesh
RQ = RandomQuadMesh
MP = MeshPlotter
PQ = PlotQuadMesh

numpts = 15
angles = range(0, 360, numpts+1)[1:end-1]
coordinates = vcat(cosd.(angles)', sind.(angles)')
boundary = [coordinates coordinates[:,1]]

hole_size = 0.5
s = hole_size/2
hole_boundary = [-s +s +s -s -s
                 -s -s +s +s -s]

# holes = [.3,.3]
# hmaxfcn = (x,y) -> 0.01 + 0.3*abs(sqrt((x-.5)^2+(y-.5)^2))

holes = [0., 0.]
p,t = RQ.polytrimesh([boundary, hole_boundary], holes,)
quad_mesh = RQ.triquad_refine(p, zeros(Int, 4, 0), t)

fig, ax = PQ.plot_mesh(quad_mesh.p, quad_mesh.t)
fig

# fig, ax = MP.plot_mesh(quad_mesh, t)
# fig