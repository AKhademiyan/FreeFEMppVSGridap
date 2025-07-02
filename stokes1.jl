using GridapMakie, GridapGmsh, Gridap
using LinearAlgebra
function main()

    # Load mesh
    model = GmshDiscreteModel("wavy.msh",has_affine_map=nothing)
    labels = get_face_labeling(model)

    # Define reference elements
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)

    # Define test spaces
    V = TestFESpace(model,reffeᵤ,labels=labels,
                    dirichlet_tags=["wall_up","wall_down"],  # Only no-slip on walls
                    conformity=:H1)
    Q = TestFESpace(model,reffeₚ,conformity=:L2)  # Remove zeromean constraint
    Y = MultiFieldFESpace([V,Q])

    # Define trial spaces with boundary conditions
    u_wall = VectorValue(0.0,0.0)  # No-slip on walls

    U = TrialFESpace(V,[u_wall,u_wall])  # Wall BCs for velocity
    P = TrialFESpace(Q)
    X = MultiFieldFESpace([U,P])

    # Define triangulation and measure
    degree = 2
    Ωₕ = Triangulation(model)
    dΩ = Measure(Ωₕ,degree)

    # Define boundary triangulation and measure
    Γ_in = BoundaryTriangulation(model, tags="inlet")
    Γ_out = BoundaryTriangulation(model, tags="outlet")
    dΓ_in = Measure(Γ_in, degree)
    dΓ_out = Measure(Γ_out, degree)

    # Define weak form
    μ = 1.0    # viscosity
    f = VectorValue(0.0,0.0)

    # Define pressure boundary conditions - Stronger pressure difference
    p_in = 100.0   # Inlet pressure
    p_out = 0.0   # Outlet pressure

    # Standard Stokes weak form with pressure boundary terms
    a((u,p),(v,q)) = ∫( μ*(∇(v)⊙∇(u)) - (∇⋅v)*p + q*(∇⋅u) )dΩ

    # Apply pressure boundary conditions directly
    l((v,q)) = ∫( v⋅f )dΩ + ∫( -p_in*(v⋅VectorValue(1.0,0.0)) )dΓ_in - ∫( -p_out*(v⋅VectorValue(1.0,0.0)) )dΓ_out

    # Create and solve the problem
    op = AffineFEOperator(a,l,X,Y)
    solve(op)
    # Solve the system
    #uh, ph = solve(op)
end

using BenchmarkTools
t = @benchmark main() samples=1 evals=1
println("Time: ", minimum(t).time/1e6, " ms")
println("Memory: ", minimum(t).memory/1024, " KiB")