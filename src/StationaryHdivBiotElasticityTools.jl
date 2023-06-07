
  const Id = TensorValue(1,0,0,1)

  default_u_ex(x)  = VectorValue(sin(π*(x[1]+x[2])), cos(π*(x[1]^2+x[2]^2)))
  default_p_ex(x)  = sin(π*x[1]+x[2])*sin(π*x[2])
  default_φE_ex(x) = sin(2π*x[1]+x[2])*sin(4π*x[2])
  default_φP_ex(x) = sin(8π*x[1]+x[2])*sin(10π*x[2])

  avXuX(v,u,μ_X,dΩ_X) = ∫(2.0*ε(v)⊙(μ_X*ε(u)))dΩ_X

  # Given the physical parameters as described in our paper we define a new
  # Biot Elast problem s.t.
  # u_E^{new} = u_E^{old}
  # u_P^{new} = u_P^{old}
  # φ_E^{new} = 1/μ_max * φ_E^{old}
  # φ_P^{new} = 1/μ_max * φ_P^{old}
  # p^{new} = 1/μ_max * p^{old}
  function scale_parameters(μ_E, μ_P, λ_E, λ_P, η, c_0, κ, α, dt; scale=true, symmetric=true)
    if (scale)
      μ_max=max(μ_E,μ_P)
    else 
      μ_max=1.0
    end
    β_E = μ_E/μ_max
    β_P = μ_P/μ_max
    ξ_E = μ_max/λ_E
    ξ_P = μ_max/λ_P
    θ1   = (α/λ_P)*μ_max
    if (symmetric)
      θ2=θ1
      δ=(c_0 + α^2/λ_P)*μ_max
      Κ   = ((dt*κ)/η)*μ_max
    else 
      θ2=(α/λ_P)
      δ=(c_0 + α^2/λ_P)
      Κ=((dt*κ)/η)
    end
    β_E, β_P, ξ_E, ξ_P, δ, θ1, θ2, Κ
  end
  # Problem definition:
  # -div(2*β_E*ε(u_E)) + div(φ_E*Id) = 0, β_E = μ_E/μ_max
  # -div(2*β_P*ε(u_P)) + div(φ_P*Id) = 0, β_P = μ_P/μ_max
  # ξ_E*φ_E + div(u_E) = 0, ξ_E = μ_max/λ_E
  # if symmetric 
     # ξ_P*φ_P - θ1*p + div(u_P) = 0, ξ_P = μ_max/λ_P
     # δ*p_P - θ2*φ_P -div(Κ*∇(p_P)) = 0, δ = (c_0 + α^2/λ_P)*μ_max, θ2 = (α/λ_P)*μ_max, Κ = ((dt*κ)/η)*μ_max
  # else 
     # ξ_P*φ_P - θ1*p + div(u_P) = 0, ξ_P = μ_max/λ_P
     # δ*p_P - θ2*φ_P -div(Κ*∇(p_P)) = 0, δ = (c_0 + α^2/λ_P), θ2 = (α/λ_P), Κ = ((dt*κ)/η)
  #    
  # With the scaling:
  # φ_E^{new} = 1/μ_max * φ_E^{old}
  # φ_P^{new} = 1/μ_max * φ_P^{old}
  # p^{new} = 1/μ_max * p^{old}
  function build_analytical_solution_terms(u_ex,p_ex,φE_ex,φP_ex,β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,ξ)
    σE_ex(x) = 2*β_E*ε(u_ex)(x)- φE_ex(x)*Id
    σP_ex(x) = 2*β_P*ε(u_ex)(x)- φP_ex(x)*Id

    σpP(x) = Κ*∇(p_ex)(x)

    buP(x) = -(∇⋅σP_ex)(x)
    buE(x) = -(∇⋅σE_ex)(x)
    bpP(x) = δ*p_ex(x) - ξ*θ2*φP_ex(x) - (∇⋅σpP)(x)

    σE_ex, σP_ex, σpP, buP, buE,  bpP
  end

 
  # ξ: can only be 1.0 or 0.0 (1: coupled pressures, 0: decoupled pressures)
  function solve_stationary_hdiv_biot_elasticity(model,
    order,
    β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,ξ;
    verbose=false,
    rtol=1.0e-06,
    fe=:raviart_thomas,
    u_ex=default_u_ex,
    p_ex=default_p_ex,
    φE_ex=default_φE_ex,
    φP_ex=default_φP_ex,
    β_U=(2.5*10^(2.0*(order-1)+1.0)),
    dt=1.0,
    solver=:lu,
    itmax=500,
    Etag="elast",
    Ptag="poroelast",
    wallEDtag="wallED",
    wallPDtag="wallPD",
    wallENtag="wallEN",
    wallPNtag="wallPN",
    dirichlet_tags_Q=[wallPDtag,wallPNtag])

    op = _assemble_stationary_hdiv_biot_elasticity(model,
                           order,
                           β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,ξ;
                           verbose=verbose,
                           fe=fe,
                           u_ex=u_ex,
                           p_ex=p_ex,
                           φE_ex=φE_ex,
                           φP_ex=φP_ex,
                           β_U=β_U,
                           dt=dt,
                           Etag=Etag,
                           Ptag=Ptag,
                           wallEDtag=wallEDtag,
                           wallPDtag=wallPDtag,
                           wallENtag=wallENtag,
                           wallPNtag=wallPNtag,
                           dirichlet_tags_Q=dirichlet_tags_Q)

    niters, uh, ph, phih = _solve_stationary_hdiv_biot_elasticity(model,
                              order,
                              op,
                              β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ;
                              verbose=verbose,
                              rtol=rtol,
                              β_U=β_U,
                              dt=dt,
                              solver=solver,
                              itmax=itmax,
                              Etag=Etag,
                              Ptag=Ptag,
                              wallEDtag=wallEDtag,
                              wallPDtag=wallPDtag)

    error_u, error_p, error_φ = _compute_errors(model,
                                               order,
                                               u_ex,
                                               p_ex,
                                               φE_ex,
                                               φP_ex,
                                               uh,ph,phih,
                                               β_U,β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,
                                               Etag=Etag,
                                               Ptag=Ptag)

    niters, error_u, error_p, error_φ, num_free_dofs(op.trial)
  end

  function _assemble_stationary_hdiv_biot_elasticity(model,order,β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,ξ;
    verbose=false,
    fe=:raviart_thomas,
    u_ex=u_ex,
    p_ex=p_ex,
    φE_ex=φE_ex,
    φP_ex=φP_ex,
    β_U=(2.5*10^(2.0*(order-1)+1.0)),
    dt=1.0,
    Etag="elast",
    Ptag="poroelast",
    wallEDtag="wallED",
    wallPDtag="wallPD",
    wallENtag="wallEN",
    wallPNtag="wallPN",
    dirichlet_tags_Q=[wallPDtag,wallPNtag])

    @assert fe == :raviart_thomas || fe==:bdm

    if (verbose)
        println("###############")
        println("order=$(order)")
        println("fe=$(fe)")
        println("β_U=$(β_U)")
        println("β_E=$(β_E)")
        println("ξ_E=$(ξ_E)")
        println("β_P=$(β_P)")
        println("ξ_P=$(ξ_P)")
        println("θ1=$(θ1)")
        println("θ2=$(θ2)")
        println("δ=$(δ)")
        println("Κ=$(Κ)")
        println("###############")
    end

    β_0 = compute_harmonic_μ_0(β_E,β_P)
    function βf(tag)
      (tag == elast) ? β_E : β_P
    end
    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))

    elast = Gridap.Geometry.get_tag_from_name(labels,Etag)
    poroelast = Gridap.Geometry.get_tag_from_name(labels,Ptag)

    Ω   = Interior(model)
    βcf=CellField(map(βf,tags),Ω)

    # Method of manufactured solutions (RHS)
    σE_ex, σP_ex, σpP, buP, buE,  bpP =
      build_analytical_solution_terms(u_ex,p_ex,φE_ex, φP_ex,β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ,ξ)

    # Reference FEs
    if (fe==:raviart_thomas)
     reffe_u   = ReferenceFE(raviart_thomas,Float64,order)
     reffe_phi = ReferenceFE(lagrangian,Float64,order)
    elseif (fe==:bdm)
     reffe_u   = ReferenceFE(bdm,Float64,order)
     reffe_phi = ReferenceFE(lagrangian,Float64,order-1)
    end
    reffe_p   = ReferenceFE(lagrangian,Float64,order)

    Ω   = Interior(model)
    Ω_E = Interior(model,tags=Etag)
    Ω_P = Interior(model,tags=Ptag)

    # Boundary triangulations and outer unit normals
    Γ_ED  = BoundaryTriangulation(model,tags=wallEDtag)
    n_ΓED = get_normal_vector(Γ_ED)
    Γ_EN  = BoundaryTriangulation(model,tags=wallENtag)
    n_ΓEN = get_normal_vector(Γ_EN)

    Γ_PD  = BoundaryTriangulation(model,tags=wallPDtag)
    n_ΓPD = get_normal_vector(Γ_PD)
    Γ_PN  = BoundaryTriangulation(model,tags=wallPNtag)
    n_ΓPN = get_normal_vector(Γ_PN)

    # Triangulation and normal on interface
    # (pointing outwards with respect to Ω_P)
    Σ_PE = InterfaceTriangulation(Ω_P,Ω_E)
    nΣ_PE  = get_normal_vector(Σ_PE)

    # Numerical integration
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΩ_P = Measure(Ω_P,degree)
    dΩ_E = Measure(Ω_E,degree)

    bdegree = 2*order
    dΓ_ED = Measure(Γ_ED,bdegree)
    dΓ_PD = Measure(Γ_PD,bdegree)
    dΓ_EN = Measure(Γ_EN,bdegree)
    dΓ_PN = Measure(Γ_PN,bdegree)


    idegree = 2*order
    dΣ_PE = Measure(Σ_PE, idegree)

    # FE spaces
    V  = TestFESpace(Ω  , reffe_u  , conformity=:HDiv, dirichlet_tags=[wallEDtag,wallPDtag])
    Q  = TestFESpace(Ω_P, reffe_p  , conformity=:H1  , dirichlet_tags=dirichlet_tags_Q)
    Z = TestFESpace(Ω, reffe_phi, conformity=:L2)


    U = TrialFESpace(V,[u_ex, u_ex])
    P = TrialFESpace(Q,p_ex)
    W = TrialFESpace(Z)
    Y = MultiFieldFESpace([V,Q,Z])
    X = MultiFieldFESpace([U,P,W])

    # Galerkin terms (no stabilisation)
    avPuP(v,u) = avXuX(v,u,β_P,dΩ_P)
    avPφP(v,φ) = ∫(-(∇⋅v)*φ)dΩ_P

    aqPpP(q,p) = ∫(δ*(p*q)+ Κ*(∇(p)⋅∇(q)))dΩ_P
    aqPφP(q,φ) = ∫((-1.0)*ξ*θ2*(φ*q))dΩ_P

    aψPφP(ψ,φ) = ∫(ξ_P*φ*ψ)dΩ_P
    aψPpP(ψ,p) = ∫(-θ1*ξ*(ψ*p))dΩ_P
    aψPuP(ψ,u) = ∫(ψ*(∇⋅u))dΩ_P

    avEuE(v,u) = avXuX(v,u,β_E,dΩ_E)
    avEφE(v,φ) = ∫(-(∇⋅v)*φ)dΩ_E

    aψEφE(ψ,φ) = ∫(ξ_E*φ*ψ)dΩ_E
    aψEuE(ψ,u) = ∫(ψ*(∇⋅u))dΩ_E

    lqP(q) = ∫(q*bpP)dΩ_P + ∫(q.⁺*((nΣ_PE.⁺)⋅σpP))dΣ_PE

    # Stabilisation terms (tangent component)
    Λ     = SkeletonTriangulation(model)
    n_Λ   = get_normal_vector(Λ)
    dΛ    = Measure(Λ   , idegree)

    Λ    = SkeletonTriangulation(model)
    Λ_E  = SkeletonTriangulation(Ω_E)
    dΛ_E = Measure(Λ_E   , idegree)
    Λ_P  = SkeletonTriangulation(Ω_P)
    dΛ_P = Measure(Λ_P   , idegree)

    n_ΛE = get_normal_vector(Λ_E)
    n_ΛP = get_normal_vector(Λ_P)

    h_e_P   = CellField(get_array(∫(1)dΛ_P), Λ_P)
    h_e_E   = CellField(get_array(∫(1)dΛ_E), Λ_E)
    h_e_PE  = CellField(get_array(∫(1)dΣ_PE), Σ_PE)
    h_e_Γ_E = CellField(get_array(∫(1)dΓ_ED), Γ_ED)
    h_e_Γ_P = CellField(get_array(∫(1)dΓ_PD), Γ_PD)

    σ(u) = 2*βcf*ε(u)
    function a_l()
      avu(v,u) = ∫(ε(v)⊙(σ(u)))dΩ - # Galerkin (bulk)
      ∫(jump(v⊗n_ΛP)⊙(mean(σ(u))))dΛ_P -
      ∫(jump(v⊗n_ΛE)⊙(mean(σ(u))))dΛ_E -
      ∫(jump(v⊗nΣ_PE)⊙(hmean(σ(u),βcf)))dΣ_PE -
      ∫((v⊗n_ΓED)⊙(σ(u)))dΓ_ED -
      ∫((v⊗n_ΓPD)⊙(σ(u)))dΓ_PD - # Galerkin (boundaries)
      ∫(jump(u⊗n_ΛP)⊙(mean(σ(v))))dΛ_P -
      ∫(jump(u⊗n_ΛE)⊙(mean(σ(v))))dΛ_E -
      ∫(jump(u⊗nΣ_PE)⊙(hmean(σ(v),βcf)))dΣ_PE -
      ∫((u⊗n_ΓED)⊙(σ(v)))dΓ_ED -
      ∫((u⊗n_ΓPD)⊙(σ(v)))dΓ_PD + # Adjoint consistency
      ∫(jump(v)⋅(2.0*β_U*β_P/h_e_P*jump(u)))dΛ_P +
      ∫(jump(v)⋅(2.0*β_U*β_E/h_e_E*jump(u)))dΛ_E +
      ∫(jump(v)⋅(2.0*β_U*β_0/h_e_PE*jump(u)))dΣ_PE +
      ∫(v⋅(2.0*β_U*β_P/h_e_Γ_P*u))dΓ_PD +
      ∫(v⋅(2.0*β_U*β_E/h_e_Γ_E*u))dΓ_ED # Stabilisation


      σuE_ex(x) = 2*β_E*ε(u_ex)(x)
      σuP_ex(x) = 2*β_P*ε(u_ex)(x)
      buuP(x) = -(∇⋅σuP_ex)(x)
      buuE(x) = -(∇⋅σuE_ex)(x)

      l(v) = ∫(v⋅buuP)dΩ_P + ∫(v⋅buuE)dΩ_E + # Galerkin (bulk)
      ∫(hmean_swap(v,βcf)⊙(jump(nΣ_PE⋅σ(u_ex))))dΣ_PE - # Galerkin boundaries
      ∫((σ(v)⋅n_ΓPD)⋅u_ex)dΓ_PD -
      ∫((σ(v)⋅n_ΓED)⋅u_ex)dΓ_ED + # Adjoint consistency
      ∫(v⋅(2.0*β_U*β_P/h_e_Γ_P*u_ex))dΓ_PD +
      ∫(v⋅(2.0*β_U*β_E/h_e_Γ_E*u_ex))dΓ_ED +
      ∫(v⋅(n_ΓEN⋅σuE_ex))dΓ_EN + # Neumann terms
      ∫(v⋅(n_ΓPN⋅σuP_ex))dΓ_PN

      return (avu, l)
    end

    avu,l=a_l()

    avφ(v,φ) = avPφP(v,φ) + avEφE(v,φ)

    φ_ex_cf=CellField(φE_ex,Ω_E)+CellField(φP_ex,Ω_P)

    ∇φP_ex(x) = ∇(φP_ex)(x)
    ∇φE_ex(x) = ∇(φE_ex)(x)
    l_φgv(v)  = ∫(∇φP_ex ⋅ v)dΩ_P +
                ∫(∇φE_ex ⋅ v)dΩ_E -
                ∫(jump(φ_ex_cf*nΣ_PE)⋅mean(v))dΣ_PE -
                ∫(jump(v⋅nΣ_PE)⋅mean(φ_ex_cf))dΣ_PE - # should be zero anyway
                ∫((v⋅n_ΓEN)*φE_ex)dΓ_EN - # Neumann terms
                ∫((v⋅n_ΓPN)*φP_ex)dΓ_PN

    l_φnew(ψ) = ∫((-ψ*θ1*ξ)*p_ex)dΩ_P +
                ∫((ξ_P)*(φP_ex*ψ))dΩ_P +
                ∫((ξ_E)*(φE_ex*ψ))dΩ_E +
                ∫(ψ*(∇⋅u_ex))dΩ_P +
                ∫(ψ*(∇⋅u_ex))dΩ_E

    # Forms
    lhs((u,p,φ),(v,q,ψ)) = avu(v,u)+avφ(v,φ)-
                           aqPpP(q,p)-aqPφP(q,φ)-
                           aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)-
                           aψEφE(ψ,φ)-aψEuE(ψ,u)

    rhs((v,q,ψ)) = l(v) + l_φgv(v) - lqP(q) - l_φnew(ψ)

    # Build affine FE operator
    op = AffineFEOperator(lhs,rhs,X,Y)
  end

  function _solve_stationary_hdiv_biot_elasticity(model,
                                     order,
                                     op,
                                     β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ;
                                     verbose=false,
                                     rtol=1.0e-06,
                                     β_U=(2.5*10^(2.0*(order-1)+1.0)),
                                     dt=1.0,
                                     solver=:lu,
                                     itmax=500,
                                     Etag="elast",
                                     Ptag="poroelast",
                                     wallEDtag="wallED",
                                     wallPDtag="wallPD")

    @assert solver==:lu || solver==:pminres

    if (solver==:pminres)
      P=_build_stationary_hdiv_biot_elasticity_riesz(model,
      order,
      op,
      β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ;
      β_U=β_U,
      dt=1.0,
      Etag=Etag,
      Ptag=Ptag,
      wallEDtag=wallEDtag,
      wallPDtag=wallPDtag)

      A = op.op.matrix
      b = op.op.vector
      X = op.trial

      r = copy(b)
      function block_stopping_condition(solver::KrylovSolver, A, b, r, tol)
        mul!(r, A, solver.x)
        r .-= b                       # r := b - Ax

        bool=true
        current=1
        for (i,X) in enumerate(op.trial)
          rv=view(r,current:current+num_free_dofs(X)-1)
          bv=view(b,current:current+num_free_dofs(X)-1)
          bool = (norm(rv) ≤ tol*norm(bv))
          @printf("[%d] ||b-Ax||/||b||: %16.7e\n",i,norm(rv)/norm(bv))
          current = current+num_free_dofs(X)
        end
        return (norm(r) ≤ tol*norm(b))
      end

      minres_callback(solver) = block_stopping_condition(solver, A, b, r, rtol)

      (xdofs, stats) = minres(A, b, M=P; itmax=itmax, verbose=1,
                              callback=minres_callback,
                              atol=0.0,
                              rtol=0.0,
                              etol=0.0)

      # Solve
      uh, ph, phih = FEFunction(X, xdofs) # solve(op)

      if (verbose)
        println("#iters=$(stats.niter)")
        println(stats)
      end

      stats.niter, uh, ph, phih
    else
      uh, ph, phih = solve(op)
      0, uh, ph, phih
    end
  end

 function _build_stationary_hdiv_biot_elasticity_riesz(model,
                                     order,
                                     op,
                                     β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ;
                                     β_U=(2.5*10^(2.0*(order-1)+1.0)),
                                     dt=1.0,
                                     Etag="elast",
                                     Ptag="poroelast",
                                     wallEDtag="wallED",
                                     wallPDtag="wallPD")
    Ω = Interior(model)
    Ω_E = Interior(model,tags=Etag)
    Ω_P = Interior(model,tags=Ptag)

    # Boundary triangulations and outer unit normals
    Γ_ED  = BoundaryTriangulation(model,tags=wallEDtag)
    n_ΓED = get_normal_vector(Γ_ED)
    Γ_PD  = BoundaryTriangulation(model,tags=wallPDtag)
    n_ΓPD = get_normal_vector(Γ_PD)

    # Numerical integration
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΩ_P = Measure(Ω_P,degree)
    dΩ_E = Measure(Ω_E,degree)

    bdegree = 2*order
    dΓ_ED = Measure(Γ_ED,bdegree)
    dΓ_PD = Measure(Γ_PD,bdegree)

    idegree = 2*order
    Λ_E  = SkeletonTriangulation(Ω_E)
    dΛ_E = Measure(Λ_E   , idegree)
    Λ_P  = SkeletonTriangulation(Ω_P)
    dΛ_P = Measure(Λ_P   , idegree)

    # Triangulation and normal on interface
    # (pointing outwards with respect to Ω_P)
    Σ_PE = InterfaceTriangulation(Ω_P,Ω_E)
    n_Σ  = get_normal_vector(Σ_PE)
    dΣ_PE = Measure(Σ_PE, idegree)

    avPuP(v,u) = avXuX(v,u,β_P,dΩ_P)
    avEuE(v,u) = avXuX(v,u,β_E,dΩ_E)
    aqPpP(q,p) = ∫(δ*(p*q)+ Κ*(∇(p)⋅∇(q)))dΩ_P

    h_e_P   = CellField(get_array(∫(1)dΛ_P), Λ_P)
    h_e_E   = CellField(get_array(∫(1)dΛ_E), Λ_E)
    h_e_PE  = CellField(get_array(∫(1)dΣ_PE), Σ_PE)
    h_e_Γ_ED = CellField(get_array(∫(1)dΓ_ED), Γ_ED)
    h_e_Γ_PD = CellField(get_array(∫(1)dΓ_PD), Γ_PD)

    β_0 = compute_harmonic_μ_0(β_E,β_P)
    riesz_map((u,p,φ),(v,q,ψ)) = avPuP(v,u) +
                                avEuE(v,u) +
                                ∫(jump(v)⋅(2.0*β_U*β_P/h_e_P*jump(u)))dΛ_P +
                                ∫(jump(v)⋅(2.0*β_U*β_E/h_e_E*jump(u)))dΛ_E +
                                ∫(jump(v)⋅(2.0*β_U*β_0/h_e_PE*jump(u)))dΣ_PE +
                                ∫(v⋅(2.0*β_U*β_E/h_e_Γ_ED*u))dΓ_ED +
                                ∫(v⋅(2.0*β_U*β_P/h_e_Γ_PD*u))dΓ_PD +
                                aqPpP(q,p) +
                                ∫((ξ_P+1.0/(2.0*β_P))*φ*ψ)dΩ_P +
                                ∫((ξ_E+1.0/(2.0*β_E))*φ*ψ)dΩ_E

    riesz_map_matrix=assemble_matrix(riesz_map,op.trial,op.test)
    P=GridapLinearSolverPreconditioner(riesz_map_matrix)      
 end 

  function _compute_errors(model,order,u_ex,p_ex,φE_ex,φP_ex,uh,ph,phih,
    β_U,β_E,β_P,ξ_E,ξ_P,θ1,θ2,δ,Κ;
                          Etag="elast",
                          Ptag="poroelast")

    degree=4*order

    Ω   = Interior(model)
    Ω_E = Interior(model,tags=Etag)
    Ω_P = Interior(model,tags=Ptag)

    # Triangulation and normal on interface
    # (pointing outwards with respect to Ω_P)
    Σ_PE = InterfaceTriangulation(Ω_P,Ω_E)

    # Numerical integration
    dΩ = Measure(Ω,degree)
    dΩ_P = Measure(Ω_P,degree)
    dΩ_E = Measure(Ω_E,degree)
    dΣ_PE = Measure(Σ_PE, degree)

    Λ  = SkeletonTriangulation(model)
    dΛ = Measure(Λ, degree)
    Λ_E=SkeletonTriangulation(Ω_E)
    Λ_P=SkeletonTriangulation(Ω_P)
    dΛ_E = Measure(Λ_E, degree)
    dΛ_P = Measure(Λ_P, degree)

    h_e_P   = CellField(get_array(∫(1)dΛ_P), Λ_P)
    h_e_E   = CellField(get_array(∫(1)dΛ_E), Λ_E)
    h_e_PE  = CellField(get_array(∫(1)dΣ_PE), Σ_PE)

    function βf(tag)
      (tag == elast) ? β_E : β_P
    end
    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))

    elast = Gridap.Geometry.get_tag_from_name(labels,Etag)
    poroelast = Gridap.Geometry.get_tag_from_name(labels,Ptag)
    βcf=CellField(map(βf,tags),Ω)

    euh        = u_ex-uh
    eph        = p_ex-ph
    ephih_E    = φE_ex-phih
    ephih_P    = φP_ex-phih

    β_0 = compute_harmonic_μ_0(β_E,β_P)
    error_u = sqrt(sum(∫(2.0*βcf*ε(euh)⊙ε(euh))dΩ  +
                       ∫(jump(euh)⋅(2.0*β_U*β_P/h_e_P*jump(euh)))dΛ_P +
                       ∫(jump(euh)⋅(2.0*β_U*β_E/h_e_E*jump(euh)))dΛ_E +
                       ∫(jump(euh)⋅(2.0*β_U*β_0/h_e_PE*jump(euh)))dΣ_PE))
    error_p = sqrt(sum(∫(δ*eph*eph  + Κ*(∇(eph)⋅∇(eph)))dΩ_P))
    error_φ = sqrt(sum(∫((1.0/β_E)*(ephih_E)*(ephih_E))dΩ_E +
                       ∫((1.0/β_P)*(ephih_P)*(ephih_P))dΩ_P))

    error_u, error_p, error_φ
  end