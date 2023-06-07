function solve_transient_hdiv_biot_elasticity(
    model,order,
    t,T,Δt,θ,u_0,p_0,φ_0,
    μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α;
    verbose=false,
    rtol=1.0e-06,
    fe=:raviart_thomas,
    β_U=(2.5*10^(2.0*(order-1)+1.0)),
    solver=:lu,
    fem_formulation=:conforming,
    solution=:analytical,
    u_ex=x->VectorValue(x[2],x[1]),
    p_ex=x->1.0,
    output_dir="transient_hdiv_biot_elasticity_solution/")

  @assert fe == :raviart_thomas || fe == :bdm
  @assert solver == :lu || solver == :pminres
  @assert fem_formulation == :conforming ||
            fem_formulation == :non_conforming
  @assert solution == :analytical || solution == :application

  rm(output_dir;force=true,recursive=true)
  mkdir(output_dir)

  if (verbose)
      println("###############")
      println("order=$(order)")
      println("fe=$(fe)")
      println("β_U=$(β_U)")
      println("μ_E=$(μ_E)")
      println("λ_E=$(λ_E)")
      println("μ_P=$(μ_P)")
      println("λ_P=$(λ_P)")
      println("η=$(η)")
      println("c_0=$(c_0)")
      println("κ=$(κ)")
      println("α=$(α)")
      println("###############")
  end

  μ_0 = compute_harmonic_μ_0(μ_E,μ_P)
  function μf(tag)
    (tag == elast) ? μ_E : μ_P
  end
  labels = get_face_labeling(model)
  tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))

  elast = Gridap.Geometry.get_tag_from_name(labels,"vessel")
  poroelast = Gridap.Geometry.get_tag_from_name(labels,"interstitium")


  if (fem_formulation == :non_conforming)
    # Reference FEs
    if (fe==:raviart_thomas)
      reffe_u = ReferenceFE(raviart_thomas,Float64,order)
      reffe_phi = ReferenceFE(lagrangian,Float64,order)
    elseif (fe==:bdm)
      reffe_u = ReferenceFE(bdm,Float64,order)
      reffe_phi = ReferenceFE(lagrangian,Float64,order-1)
    end
    reffe_p = ReferenceFE(lagrangian,Float64,order)
  else
    # Reference FEs
    reffe_u   = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffe_p   = ReferenceFE(lagrangian,Float64,order-1)
    reffe_phi = ReferenceFE(lagrangian,Float64,order-2)
  end

  Ω   = Interior(model)
  Ω_E = Interior(model,tags="vessel")
  Ω_P = Interior(model,tags="interstitium")

  μcf=CellField(map(μf,tags),Ω)
  
  # Boundary triangulations and outer unit normals

  if (solution==:analytical)
    # Boundary triangulations and outer unit normals
    Γ_E = BoundaryTriangulation(model,tags=["free_vessel","clamped_vessel", "traction"])
    Γ_P = BoundaryTriangulation(model,tags=["free_interstitium", "clamped_interstitium"])

    Γ_PN = BoundaryTriangulation(model,tags=["free_interstitium"])
    Γ_PD = BoundaryTriangulation(model,tags=["clamped_interstitium"])
    n_ΓPN = get_normal_vector(Γ_PN)
    n_ΓPD = get_normal_vector(Γ_PD)

    Γ_EN = BoundaryTriangulation(model,tags=["free_vessel","traction"])
    Γ_ED = BoundaryTriangulation(model,tags=["clamped_vessel"])
    n_ΓEN = get_normal_vector(Γ_EN)
    n_ΓED = get_normal_vector(Γ_ED)
  else
    Γ_E = BoundaryTriangulation(model,tags=["free_vessel","clamped_vessel", "traction"])
    Γ_E_traction = BoundaryTriangulation(model,tags=["traction"])
    Γ_EN = BoundaryTriangulation(model,tags=["free_vessel","traction"])
    Γ_ED = BoundaryTriangulation(model,tags=["clamped_vessel"])
    n_ΓE = get_normal_vector(Γ_E)
    n_ΓE_traction = get_normal_vector(Γ_E_traction)
    n_ΓEN = get_normal_vector(Γ_EN)
    n_ΓED = get_normal_vector(Γ_ED)

    Γ_P = BoundaryTriangulation(model,tags=["free_interstitium", "clamped_interstitium"])
    Γ_PN = BoundaryTriangulation(model,tags=["free_interstitium"])
    Γ_PD = BoundaryTriangulation(model,tags=["clamped_interstitium"])
    n_ΓP = get_normal_vector(Γ_P)
    n_ΓPN = get_normal_vector(Γ_PN)
    n_ΓPD = get_normal_vector(Γ_PD)
  end

  # Triangulation and normal on interface
  # (pointing outwards with respect to Ω_P)
  Σ_PE = InterfaceTriangulation(Ω_P,Ω_E)
  n_Σ  = get_normal_vector(Σ_PE)

  # Numerical integration
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΩ_P = Measure(Ω_P,degree)
  dΩ_E = Measure(Ω_E,degree)

  bdegree = 2*order
  if (solution==:analytical)
    dΓ_E = Measure(Γ_E,bdegree)
    dΓ_P = Measure(Γ_P,bdegree)
    dΓ_ED = Measure(Γ_ED,bdegree)
    dΓ_EN = Measure(Γ_EN,bdegree)
    dΓ_PD = Measure(Γ_PD,bdegree)
    dΓ_PN = Measure(Γ_PN,bdegree)
  else
    dΓ_E = Measure(Γ_E,bdegree)
    dΓ_P = Measure(Γ_P,bdegree)
    dΓ_ED = Measure(Γ_ED,bdegree)
    dΓ_EN = Measure(Γ_EN,bdegree)
    dΓ_E_traction = Measure(Γ_E_traction,bdegree)
    dΓ_PD = Measure(Γ_PD,bdegree)
    dΓ_PN = Measure(Γ_PN,bdegree)
  end

  idegree = 2*order
  dΣ_PE = Measure(Σ_PE, idegree)

  # Stabilisation terms (tangent component)
  Λ   = SkeletonTriangulation(model)
  Λ_E = SkeletonTriangulation(Ω_E)
  Λ_P = SkeletonTriangulation(Ω_P)

  n_ΛE = get_normal_vector(Λ_E)
  n_ΛP = get_normal_vector(Λ_P)
  nΣ_PE  = get_normal_vector(Σ_PE)

  dΛ_E  = Measure(Λ_E , idegree)
  dΛ_P  = Measure(Λ_P , idegree)
  dΛ    = Measure(Λ   , idegree)

  h_e_P   = CellField(get_array(∫(1)dΛ_P), Λ_P)
  h_e_E   = CellField(get_array(∫(1)dΛ_E), Λ_E)
  h_e_PE  = CellField(get_array(∫(1)dΣ_PE), Σ_PE)

  if (solution==:analytical)
    h_e_Γ_ED = CellField(get_array(∫(1)dΓ_ED), Γ_ED)
    h_e_Γ_PD = CellField(get_array(∫(1)dΓ_PD), Γ_PD)
    h_e_Γ_P = CellField(get_array(∫(1)dΓ_P), Γ_P)
    h_e_Γ_E = CellField(get_array(∫(1)dΓ_E), Γ_E)

  else
    h_e_Γ_ED = CellField(get_array(∫(1)dΓ_ED), Γ_ED)
    h_e_Γ_PD = CellField(get_array(∫(1)dΓ_PD), Γ_PD)
    h_e_Γ_P = CellField(get_array(∫(1)dΓ_P), Γ_P)
    h_e_Γ_E = CellField(get_array(∫(1)dΓ_E), Γ_E)
  end

  u_ex_t=nothing
  p_ex_t=nothing
  if (solution==:analytical)
    φE_ex(x) = -λ_E*(∇⋅u_ex)(x)
    φE_ex_t(x,t::Real) = φE_ex(x)
    φE_ex_t(t::Real) = (x) -> φE_ex_t(x,t)

    φP_ex(x) = α*p_ex(x)-λ_P*(∇⋅u_ex)(x)
    φP_ex_t(x,t::Real) = φP_ex(x)
    φP_ex_t(t::Real) = (x) -> φP_ex_t(x,t)

    σE_ex(x) = 2*μ_E*ε(u_ex)(x)- φE_ex(x)*Id
    σE_ex_t(x,t::Real) = σE_ex(x)
    σE_ex_t(t::Real) = (x) -> σE_ex_t(x,t)
    σP_ex(x) = 2*μ_P*ε(u_ex)(x)- φP_ex(x)*Id
    σP_ex_t(x,t::Real) = σP_ex(x)
    σP_ex_t(t::Real) = (x) -> σP_ex_t(x,t)

    σpP(x) = κ/η*∇(p_ex)(x)
    σpP_t(x,t::Real) = σpP(x)
    σpP_t(t::Real) = (x) -> σpP_t(x,t)

    buP(x) = -(∇⋅σP_ex)(x)
    buP_t(x,t::Real) = buP(x)
    buP_t(t::Real) = (x) -> buP_t(x,t)

    buE(x) = -(∇⋅σE_ex)(x)
    buE_t(x,t::Real) = buE(x)
    buE_t(t::Real) = (x) -> buE_t(x,t)

    σuE_ex(x) = 2*μ_E*ε(u_ex)(x)
    σuE_ex_t(x,t::Real) = σuE_ex(x)
    σuE_ex_t(t::Real) = (x) -> σuE_ex_t(x,t)

    σuP_ex(x) = 2*μ_P*ε(u_ex)(x)
    σuP_ex_t(x,t::Real) = σuP_ex(x)
    σuP_ex_t(t::Real) = (x) -> σuP_ex_t(x,t)

    buuP(x) = -(∇⋅σuP_ex)(x)
    buuP_t(x,t::Real) = buuP(x)
    buuP_t(t::Real) = (x) -> buuP_t(x,t)

    buuE(x) = -(∇⋅σuE_ex)(x)
    buuE_t(x,t::Real) = buuE(x)
    buuE_t(t::Real) = (x) -> buuE_t(x,t)

    ∇φP_ex(x) = ∇(φP_ex)(x)
    ∇φP_ex(x,t::Real)=∇φP_ex(x)
    ∇φP_ex_t(t::Real)=(x) -> ∇φP_ex(x,t)

    ∇φE_ex(x) = ∇(φE_ex)(x)
    ∇φE_ex(x,t::Real)=∇φE_ex(x)
    ∇φE_ex_t(t::Real)=(x) -> ∇φE_ex(x,t)

    u_ex_t(x,t::Real) = u_ex(x)
    u_ex_t(t::Real)   = (x) -> u_ex_t(x,t)
    p_ex_t(x,t::Real) = p_ex(x)
    p_ex_t(t::Real)   = (x) -> p_ex_t(x,t)

    bpP_t(t::Real) = (x) -> (c_0 + α^2/λ_P)*∂t(p_ex_t)(t)(x) - α/λ_P*∂t(φP_ex_t)(t)(x) - (∇⋅(σpP_t(t)))(x)
  end

  u_conformity=:H1
  if (fem_formulation == :non_conforming)
    u_conformity=:HDiv
  end
  if (solution == :analytical)
    V  = TestFESpace(Ω, reffe_u, conformity=u_conformity, dirichlet_tags=["clamped_interstitium", "clamped_vessel"])
    Q  = TestFESpace(Ω_P, reffe_p, conformity=:H1, dirichlet_tags=["clamped_interstitium"])
    Z  = TestFESpace(Ω, reffe_phi, conformity=:L2)
    U  = TransientTrialFESpace(V, [u_ex_t, u_ex_t])
    P  = TransientTrialFESpace(Q, [p_ex_t])
    W  = TrialFESpace(Z)
    Y = MultiFieldFESpace([V,Q,Z])
    X = TransientMultiFieldFESpace([U,P,W])
  else # non-analytical solution (application)
    V  = TestFESpace(Ω, reffe_u, conformity=u_conformity, dirichlet_tags=["clamped_interstitium", "clamped_vessel"])
    Q  = TestFESpace(Ω_P, reffe_p, conformity=:H1, dirichlet_tags="clamped_interstitium")
    Z  = TestFESpace(Ω, reffe_phi, conformity=:L2)
    # Dirichlet BC
    u_diri(x) = VectorValue(0,0)
    #u_diri(t::Real) = x -> u_diri(x,t)
    U  = TrialFESpace(V, [u_diri, u_diri])
    W  = TrialFESpace(Z)
    # g1(x,t::Real) = 1.0e4*sin(2*π*t)
    # g1(t::Real)   = x -> g1(x,t)

    p_diri(x,t::Real) = 0
    p_diri(t::Real) = x ->  p_diri(x,t)
    P = TransientTrialFESpace(Q, p_diri)
    Y = MultiFieldFESpace([V,Q,Z])
    X = TransientMultiFieldFESpace([U,P,W])
  end

  # Defining forms common either all or some of the methods

  # Coefficient to sym the operator
  Δtθ = Δt * θ

  # Galerkin terms (no stabilisation)
  avPuP(v,u) = ∫(2.0*ε(v)⊙(μ_P*ε(u)))dΩ_P
  avPφP(v,φ) = ∫(-(∇⋅v)*φ)dΩ_P

  # Define the two terms with time derivatives
  aqPpPt(q,dp) = ∫(Δtθ*(-c_0-α^2/λ_P)*(dp*q))dΩ_P
  aqPφPt(q,dφ) = ∫(-Δtθ*α/λ_P*(dφ*q))dΩ_P

  aqPpP(q,p) = ∫(Δtθ*κ/η*(∇(p)⋅∇(q)))dΩ_P
  aψPφP(ψ,φ) = ∫(1.0/λ_P*φ*ψ)dΩ_P
  aψPpP(ψ,p) = ∫(-α/λ_P*(ψ*p))dΩ_P
  aψPuP(ψ,u) = ∫(ψ*(∇⋅u))dΩ_P

  avEuE(v,u) = ∫(2.0*ε(v)⊙(μ_E*ε(u)))dΩ_E
  avEφE(v,φ) = ∫(-(∇⋅v)*φ)dΩ_E

  aψEφE(ψ,φ) = ∫(1.0/λ_E*φ*ψ)dΩ_E
  aψEuE(ψ,u) = ∫(ψ*(∇⋅u))dΩ_E

  # Function definition required for application problem
  # Elasticity momentum equation Neumann BC data
  g2(x,t::Real) = 1e3*sin(2*π*t)
  g2(t::Real) = x -> g2(x,t)

  function generate_application_transient_operator_forms_non_conforming()
    σ(u) = 2*μcf*ε(u)

    avu(v,u) = ∫(ε(v)⊙(σ(u)))dΩ - # Galerkin (bulk)
              ∫(jump(v⊗n_ΛP)⊙(mean(σ(u))))dΛ_P -
              ∫(jump(v⊗n_ΛE)⊙(mean(σ(u))))dΛ_E -
              ∫(jump(v⊗nΣ_PE)⊙(hmean(σ(u),μcf)))dΣ_PE -
              ∫((v⊗n_ΓED)⊙(σ(u)))dΓ_ED -
              ∫((v⊗n_ΓPD)⊙(σ(u)))dΓ_PD - # Galerkin (boundaries)
              ∫(jump(u⊗n_ΛP)⊙(mean(σ(v))))dΛ_P -
              ∫(jump(u⊗n_ΛE)⊙(mean(σ(v))))dΛ_E -
              ∫(jump(u⊗nΣ_PE)⊙(hmean(σ(v),μcf)))dΣ_PE -
              ∫((u⊗n_ΓED)⊙(σ(v)))dΓ_ED -
              ∫((u⊗n_ΓPD)⊙(σ(v)))dΓ_PD + # Adjoint consistency
              ∫(jump(v)⋅(2.0*β_U*μ_P/h_e_P*jump(u)))dΛ_P +
              ∫(jump(v)⋅(2.0*β_U*μ_E/h_e_E*jump(u)))dΛ_E +
              ∫(jump(v)⋅(2.0*β_U*μ_0/h_e_PE*jump(u)))dΣ_PE +
              ∫(v⋅(2.0*β_U*μ_P/h_e_Γ_PD*u))dΓ_PD +
              ∫(v⋅(2.0*β_U*μ_E/h_e_Γ_ED*u))dΓ_ED # Stabilisation

    # Forms
    avφ(v,φ) = avPφP(v,φ) + avEφE(v,φ)

    lhs(t,(u,p,φ),(v,q,ψ)) = avu(v,u)+avφ(v,φ)-
                            aqPpP(q,p)-
                            aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)-
                            aψEφE(ψ,φ)-aψEuE(ψ,u)

    m(t,(dtu,dtp,dtφ),(v,q,ψ)) = aqPpPt(q,dtp)-aqPφPt(q,dtφ)
    b(t,(v,q,ψ)) = ∫(g2(t)*(n_ΓE_traction⋅v))dΓ_E_traction
    m,lhs,b
  end

  function generate_analytical_transient_operator_forms_non_conforming()
    σ(u) = 2*μcf*ε(u)
    avu(v,u) = ∫(ε(v)⊙(σ(u)))dΩ - # Galerkin (bulk)
      ∫(jump(v⊗n_ΛP)⊙(mean(σ(u))))dΛ_P -
      ∫(jump(v⊗n_ΛE)⊙(mean(σ(u))))dΛ_E -
      ∫(jump(v⊗nΣ_PE)⊙(hmean(σ(u),μcf)))dΣ_PE -
      ∫((v⊗n_ΓED)⊙(σ(u)))dΓ_ED -
      ∫((v⊗n_ΓPD)⊙(σ(u)))dΓ_PD - # Galerkin (boundaries)
      ∫(jump(u⊗n_ΛP)⊙(mean(σ(v))))dΛ_P -
      ∫(jump(u⊗n_ΛE)⊙(mean(σ(v))))dΛ_E -
      ∫(jump(u⊗nΣ_PE)⊙(hmean(σ(v),μcf)))dΣ_PE -
      ∫((u⊗n_ΓED)⊙(σ(v)))dΓ_ED -
      ∫((u⊗n_ΓPD)⊙(σ(v)))dΓ_PD + # Adjoint consistency
      ∫(jump(v)⋅(2.0*β_U*μ_P/h_e_P*jump(u)))dΛ_P +
      ∫(jump(v)⋅(2.0*β_U*μ_E/h_e_E*jump(u)))dΛ_E +
      ∫(jump(v)⋅(2.0*β_U*μ_0/h_e_PE*jump(u)))dΣ_PE +
      ∫(v⋅(2.0*β_U*μ_P/h_e_Γ_PD*u))dΓ_PD +
      ∫(v⋅(2.0*β_U*μ_E/h_e_Γ_ED*u))dΓ_ED # Stabilisation


    φ_ex_cf=CellField(φE_ex,Ω_E)+CellField(φP_ex,Ω_P)
    l(t,v) = ∫(v⋅buuP_t(t))dΩ_P + ∫(v⋅buuE_t(t))dΩ_E +
      ∫(hmean_swap(v,μcf)⊙(jump(nΣ_PE⋅σ(u_ex_t(t)))))dΣ_PE -
      ∫((σ(v)⋅n_ΓPD)⋅u_ex_t(t))dΓ_PD -
      ∫((σ(v)⋅n_ΓED)⋅u_ex_t(t))dΓ_ED + # Adjoint consistency
      ∫(v⋅(2.0*β_U*μ_P/h_e_Γ_PD*u_ex_t(t)))dΓ_PD +
      ∫(v⋅(2.0*β_U*μ_E/h_e_Γ_ED*u_ex_t(t)))dΓ_ED +
      ∫(v⋅(n_ΓEN⋅σuE_ex_t(t)))dΓ_EN + # Neumann terms
      ∫(v⋅(n_ΓPN⋅σuP_ex_t(t)))dΓ_PN

    l_φgv(t,v)  = ∫(∇φP_ex_t(t) ⋅ v)dΩ_P +
      ∫(∇φE_ex_t(t) ⋅ v)dΩ_E -
      ∫(jump(φ_ex_cf*nΣ_PE)⋅hmean_swap(v,μcf))dΣ_PE -
      ∫(jump(v⋅nΣ_PE)⋅mean(φ_ex_cf))dΣ_PE - # should be zero anyway
      ∫((v⋅n_ΓEN)*φE_ex_t(t))dΓ_EN - # Neumann terms
      ∫((v⋅n_ΓPN)*φP_ex_t(t))dΓ_PN

    lqP(t,q) = ∫(q*bpP_t(t))dΩ_P + ∫(q.⁺*((n_Σ.⁺)⋅σpP_t(t)))dΣ_PE

    avφ(v,φ) = avPφP(v,φ) + avEφE(v,φ)


    lhs(t,(u,p,φ),(v,q,ψ)) = avu(v,u)+avφ(v,φ)-
                                 aqPpP(q,p)-
                                 aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)-
                                 aψEφE(ψ,φ)-aψEuE(ψ,u)

    m(t,(dtu,dtp,dtφ),(v,q,ψ)) = aqPpPt(q,dtp)-aqPφPt(q,dtφ)

    b(t,(v,q,ψ)) = l(t,v) + l_φgv(t,v) - lqP(t,q)
    m,lhs,b
  end

  function generate_application_transient_operator_forms_conforming()
    # Forms
    σ(u) = 2*μcf*ε(u)

    avu(v,u) = ∫(ε(v)⊙(σ(u)))dΩ

    # Forms
    avφ(v,φ) = avPφP(v,φ) + avEφE(v,φ)

    lhs(t,(u,p,φ),(v,q,ψ)) = avu(v,u)+avφ(v,φ)-
                            aqPpP(q,p)-
                            aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)-
                            aψEφE(ψ,φ)-aψEuE(ψ,u)

    m(t,(dtu,dtp,dtφ),(v,q,ψ)) = aqPpPt(q,dtp)-aqPφPt(q,dtφ)
    b(t,(v,q,ψ)) = ∫(g2(t)*(n_ΓE_traction⋅v))dΓ_E_traction
    m,lhs,b
  end

  function generate_analytical_transient_operator_forms_conforming()
    φE_ex(x) = -λ_E*(∇⋅u_ex)(x)
    φP_ex(x) = α*p_ex(x)-λ_P*(∇⋅u_ex)(x)
    σE_ex(x) = 2*μ_E*ε(u_ex)(x)- φE_ex(x)*Id
    σP_ex(x) = 2*μ_P*ε(u_ex)(x)- φP_ex(x)*Id
    σpP_ex(x) = κ/η*∇(p_ex)(x)
    buP_ex(x) = -(∇⋅σP_ex)(x)
    buE_ex(x) = -(∇⋅σE_ex)(x)

    u_ex_t(x,t::Real) = u_ex(x)
    u_ex_t(t::Real)   = (x) -> u_ex_t(x,t)
    p_ex_t(x,t::Real) = p_ex(x)
    p_ex_t(t::Real)   = (x) -> p_ex_t(x,t)

    φE_ex_t(x,t::Real) = φE_ex(x)
    φE_ex_t(t::Real) = (x) -> φE_ex_t(x,t)
    φP_ex_t(x,t::Real) = φP_ex(x)
    φP_ex_t(t::Real) = (x) -> φP_ex_t(x,t)
    σE_ex_t(x,t::Real) = σE_ex(x)
    σE_ex_t(t::Real) = (x) -> σE_ex_t(x,t)
    σP_ex_t(x,t::Real) = σP_ex(x)
    σP_ex_t(t::Real) = (x) -> σP_ex_t(x,t)
    σpP_ex_t(x,t::Real) = σpP_ex(x)
    σpP_ex_t(t::Real) = (x) -> σpP_ex_t(x,t)
    buP_ex_t(x,t::Real) = buP_ex(x)
    buP_ex_t(t::Real) = (x) -> buP_ex_t(x,t)
    buE_ex_t(x,t::Real) = buE_ex(x)
    buE_ex_t(t::Real) = (x) -> buE_ex_t(x,t)
    bpP_ex_t(t::Real) = (x) -> (c_0 + α^2/λ_P)*∂t(p_ex_t)(t)(x) - α/λ_P*∂t(φP_ex_t)(t)(x) - (∇⋅(σpP_ex_t(t)))(x)

    # Forms required for analytical solutions
    lvP(t,v) = ∫(v⋅buP_ex_t(t))dΩ_P + ∫(v.⁺⋅((n_Σ.⁺)⋅σP_ex_t(t)))dΣ_PE
    lvE(t,v) = ∫(v⋅buE_ex_t(t))dΩ_E + ∫(v.⁻⋅((n_Σ.⁻)⋅σE_ex_t(t)))dΣ_PE
    lqP(t,q) = ∫(q*bpP_ex_t(t))dΩ_P + ∫(q.⁺*((n_Σ.⁺)⋅σpP_ex_t(t)))dΣ_PE

    # Forms
    lhs(t,(u,p,φ),(v,q,ψ)) = avPuP(v,u)+avPφP(v,φ)-
                             aqPpP(q,p)-
                             aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)+
                             avEuE(v,u)+avEφE(v,φ)-
                             aψEφE(ψ,φ)-aψEuE(ψ,u)
    m(t,(dtu,dtp,dtφ),(v,q,ψ)) = aqPpPt(q,dtp)-aqPφPt(q,dtφ)
    b(t,(v,q,ψ)) = lvP(t,v)+lvE(t,v)-lqP(t,q) +
                  ∫(v⋅(n_ΓEN⋅σuE_ex_t(t)))dΓ_EN +
                  ∫(v⋅(n_ΓPN⋅σuP_ex_t(t)))dΓ_PN -
                  ∫((v⋅n_ΓEN)*φE_ex_t(t))dΓ_EN -
                  ∫((v⋅n_ΓPN)*φP_ex_t(t))dΓ_PN
    m,lhs,b
  end

  if (fem_formulation==:non_conforming)
    if (solution==:analytical)
      m,lhs,b=generate_analytical_transient_operator_forms_non_conforming()
    else
      m,lhs,b=generate_application_transient_operator_forms_non_conforming()
    end
  else
    if (solution==:analytical)
      m,lhs,b=generate_analytical_transient_operator_forms_conforming()
    else
      m,lhs,b=generate_application_transient_operator_forms_conforming()
    end
  end

  # Build affine FE operator
  op = TransientAffineFEOperator(m,lhs,b,X,Y)

  if (solver==:lu)
    linear_solver = LUSolver()
  else
    riesz_map=nothing
    if (solution==:application && fem_formulation==:non_conforming)
          rm_app_nc((u,p,φ),(v,q,ψ)) = avPuP(v,u) +
                                      avEuE(v,u) +
                                      ∫(jump(v)⋅(2.0*β_U*μ_P/h_e_P*jump(u)))dΛ_P +
                                      ∫(jump(v)⋅(2.0*β_U*μ_E/h_e_E*jump(u)))dΛ_E +
                                      ∫(jump(v)⋅(2.0*β_U*μ_0/h_e_PE*jump(u)))dΣ_PE +
                                      ∫(v⋅(2.0*β_U*μ_P/h_e_Γ_P*u))dΓ_PD +
                                      ∫(v⋅(2.0*β_U*μ_E/h_e_Γ_E*u))dΓ_ED +
                                      aqPpP(q,p) + ∫((c_0+α^2/λ_P)*(p*q))dΩ_P +
                                      ∫((1.0/λ_P+1.0/(2.0*μ_P))*φ*ψ)dΩ_P +
                                      ∫((1.0/λ_E+1.0/(2.0*μ_E))*φ*ψ)dΩ_E
          riesz_map_matrix=assemble_matrix(rm_app_nc,MultiFieldFESpace(TrialFESpace.([V,Q,Z])),Y)
    elseif (solution==:analytical && fem_formulation==:non_conforming)
          rm_analytical_nc((u,p,φ),(v,q,ψ)) = avPuP(v,u) + avEuE(v,u) +
                                 ∫(jump(v)⋅(2.0*β_U*μ_P/h_e_P*jump(u)))dΛ_P +
                                 ∫(jump(v)⋅(2.0*β_U*μ_E/h_e_E*jump(u)))dΛ_E +
                                 ∫(jump(v)⋅(2.0*β_U*μ_0/h_e_PE*jump(u)))dΣ_PE +
                                 ∫(v⋅(2.0*β_U*μ_P/h_e_Γ_P*u))dΓ_PD +
                                 ∫(v⋅(2.0*β_U*μ_E/h_e_Γ_E*u))dΓ_ED +
                                 aqPpP(q,p) + ∫((c_0+α^2/λ_P)*(p*q))dΩ_P +
                                 ∫((1.0/λ_P+1.0/(2.0*μ_P))*φ*ψ)dΩ_P +
                                 ∫((1.0/λ_E+1.0/(2.0*μ_E))*φ*ψ)dΩ_E
          riesz_map_matrix=assemble_matrix(rm_analytical_nc,MultiFieldFESpace(TrialFESpace.([V,Q,Z])),Y)
    elseif (solution==:application && fem_formulation==:conforming)
          rm_app_c((u,p,φ),(v,q,ψ)) = avPuP(v,u) +
                                      avEuE(v,u) +
                                      aqPpP(q,p) + ∫((1.0/Δt)*(c_0+α^2/λ_P)*(p*q))dΩ_P +
                                      ∫((1.0/λ_P+1.0/(2.0*μ_P))*φ*ψ)dΩ_P +
                                      ∫((1.0/λ_E+1.0/(2.0*μ_E))*φ*ψ)dΩ_E
          riesz_map_matrix=assemble_matrix(rm_app_c,MultiFieldFESpace(TrialFESpace.([V,Q,Z])),Y)
    elseif (solution==:analytical && fem_formulation==:conforming)
          rm_analytical_c((u,p,φ),(v,q,ψ)) = avPuP(v,u) +
                                      avEuE(v,u) +
                                      aqPpP(q,p) + ∫((1.0/Δt)*(c_0+α^2/λ_P)*(p*q))dΩ_P +
                                      ∫((1.0/λ_P+1.0/(2.0*μ_P))*φ*ψ)dΩ_P +
                                      ∫((1.0/λ_E+1.0/(2.0*μ_E))*φ*ψ)dΩ_E
          riesz_map_matrix=assemble_matrix(rm_analytical_c,MultiFieldFESpace(TrialFESpace.([V,Q,Z])),Y)

    end
    M=GridapLinearSolverPreconditioner(riesz_map_matrix)
    linear_solver=PMinResSolver(M;itmax=500,rtol=rtol,verbose=1)
  end

  ode_solver = ThetaMethod(linear_solver,Δt,θ)
  if (solution==:application)
    u_θ = interpolate_everywhere(u_0,U(t))
    p_θ = interpolate_everywhere(p_0,P(t))
    φ_θ = interpolate_everywhere(φ_0,W(t))
    x_θ = interpolate_everywhere([u_θ,p_θ,φ_θ],X(t));
  else
    φ_ex_cf=CellField(φE_ex_t(t),Ω_E)+CellField(φP_ex_t(t),Ω_P)
    u_θ = interpolate_everywhere(u_ex_t(t),U(t))
    p_θ = interpolate_everywhere(p_ex_t(t),P(t))
    φ_θ = interpolate_everywhere(φ_ex_cf,W)
    x_θ = interpolate_everywhere([u_θ,p_θ,φ_θ],X(t));
  end

  sol = solve(ode_solver,op,x_θ,t,T)

  createpvd("$(output_dir)/transient_solution_biot_elasticity") do pvd
    if (solution==:analytical)
      error_u = sqrt(sum(∫((u_ex-u_θ)⋅(u_ex-u_θ))*dΩ))
      error_p = sqrt(sum(∫((p_ex-p_θ)⋅(p_ex-p_θ))*dΩ))
      error_φ = sqrt(sum(∫((φ_ex_cf-φ_θ)⋅(φ_ex_cf-φ_θ))*dΩ))
      println("-----")
      println("error_u ($t)=$(error_u)")
      println("error_p ($t)=$(error_p)")
      println("error_φ ($t)=$(error_φ)")
      println("-----")
    end
    pvd[t] = createvtk(Ω,
                      "$(output_dir)/transient_solution_biot_elasticity_$t"*".vtu",
                       cellfields=["u_n"=>u_θ,"p_n"=>p_θ,"φ_n"=>φ_θ])
    for (x_n,t_n) in sol
      (u_n,p_n,φ_n) = x_n
      if (solution==:analytical)
         error_u = sqrt(sum(∫((u_ex-u_n)⋅(u_ex-u_n))*dΩ))
         error_p = sqrt(sum(∫((p_ex-p_n)⋅(p_ex-p_n))*dΩ))
         error_φ = sqrt(sum(∫((φ_ex_cf-φ_n)⋅(φ_ex_cf-φ_n))*dΩ))
         println("-----")
         println("error_u ($t_n)=$(error_u)")
         println("error_p ($t_n)=$(error_p)")
         println("error_φ ($t_n)=$(error_φ)")
         println("-----")
      else 
        if (verbose)
          println("T=($t_n) completed!")
        end   
      end
      pvd[t_n] = createvtk(Ω,
                          "$(output_dir)/transient_solution_biot_elasticity_$t_n"*".vtu",
                          cellfields=["u_n"=>u_n,"p_n"=>p_n,"φ_n"=>φ_n])
    end
  end
end