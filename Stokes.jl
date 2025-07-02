###############################################################################
# wavy_stokes.jl  –  put *everything* below this header exactly in this order
###############################################################################

using Gmsh: gmsh                       # ① LOADS Gmsh.jl   → creates the constant `gmsh`
println("→ Gmsh library path: ", gmsh.get_install_message())
println("→ Gmsh version: ", gmsh.version())
using Gridap, GridapGmsh
using LinearAlgebra

# ─────────────────────────── parameters ─────────────────────────────
const L = 1.0
const α = 0.2
const ε = 0.8
const nnL = 300
const nnH = 40
const order = 2
const H = α*L
const k = 2π/L
const μ = 1.0
const p_in = 1.0
const p_out = 0.0

# ────────────────────────── geometry  (GEO API) ─────────────────────
gmsh.initialize()
gmsh.model.add("wavy_channel")
gp = gmsh.model.geo                         # use the Geo kernel (works w/o OCC)

xs   = range(0, L; length = nnL+1)
hmax = max(L/nnL, 2H/nnH)

p_bot = [gp.addPoint(x, -H*(1-ε*cos(k*x)), 0, hmax) for x in xs]
p_top = [gp.addPoint(x,  H*(1-ε*cos(k*x)), 0, hmax) for x in xs]

l_bot   = gp.addSpline(p_bot)
l_right = gp.addLine(p_top[end], p_bot[end])
l_top   = gp.addSpline(p_top)
l_left  = gp.addLine(p_bot[1],  p_top[1])

loop  = gp.addCurveLoop([l_bot, l_right, l_top, l_left])
surf  = gp.addPlaneSurface([loop])
gp.synchronize()                             # ← Geo sync, not OCC

gmsh.model.addPhysicalGroup(1,[l_bot], 1);   gmsh.model.setPhysicalName(1,1,"wall_down")
gmsh.model.addPhysicalGroup(1,[l_right],2); gmsh.model.setPhysicalName(1,2,"outlet")
gmsh.model.addPhysicalGroup(1,[l_top], 3);   gmsh.model.setPhysicalName(1,3,"wall_up")
gmsh.model.addPhysicalGroup(1,[l_left], 4);  gmsh.model.setPhysicalName(1,4,"inlet")
gmsh.model.addPhysicalGroup(2,[surf], 5);    gmsh.model.setPhysicalName(2,5,"domain")

gmsh.option.setNumber("Mesh.ElementOrder", order)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", hmax)
gmsh.model.mesh.generate(2)
# ───────────────────────── Gridap model & spaces ─────────────────────
model  = GmshDiscreteModel(gmsh.model)
labels = get_face_labeling(model)

orderu = order
orderp = order-1
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},orderu)
reffeₚ = ReferenceFE(lagrangian,Float64,orderp; space=:P)

V = TestFESpace(model,reffeᵤ; labels, dirichlet_tags=["wall_up","wall_down"],
                conformity=:H1)
Q = TestFESpace(model,reffeₚ; conformity=:L2)
Y = MultiFieldFESpace([V,Q])

u_wall = VectorValue(0.0,0.0)
U = TrialFESpace(V,[u_wall,u_wall])
P = TrialFESpace(Q)
X = MultiFieldFESpace([U,P])

Ω  = Triangulation(model)
Γin  = BoundaryTriangulation(model,tags="inlet")
Γout = BoundaryTriangulation(model,tags="outlet")
dΩ    = Measure(Ω, order)
dΓin  = Measure(Γin, order)
dΓout = Measure(Γout, order)

a((u,p),(v,q)) = ∫( μ*(∇(v)⊙∇(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ
l((v,q)) = ∫(                            )dΩ +
           ∫( -(p_in)*(v⋅VectorValue(1,0)) )dΓin +
           ∫( -(p_out)*(v⋅VectorValue(1,0)) )dΓout

println("► Assembling & solving Stokes system …")
op = AffineFEOperator(a,l,X,Y)
uh,ph = solve(op)
println("   Done.")

# ─────────────────── Post-process pressure & vorticity ───────────────
min_p   = minimum(get_free_dof_values(ph))
offset  = abs(min_p) + 0.1
ph_shift = ph + offset

# Vorticity as a cell field
vort = evaluate(
  CellField((x, u) -> ∂(u[2], 1) - ∂(u[1], 2), uh),
  Ω
)

println("   Pressure shifted by ",offset," to ensure positivity.")

# ──────────────────────── Export VTK file ────────────────────────────
writevtk(Ω,"stokes_pressure_driven",
         cellfields=["velocity"=>uh,
                     "pressure"=>ph_shift,
                     "vorticity"=>vort])
println("► Results written to stokes_pressure_driven.vtu")

gmsh.finalize()
#######################################################################

