module Driver
  # Add here whichever modules/packages your
  # driver function may need
  using HdivBiotElasticityPaper
  using Gridap

  function driver(problem,solver,fe,nk,order,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,scale)
    @assert problem in (:square,:vessel_bl,:vessel_t)
    if (problem==:square) 
      model=generate_model_unit_square_biot_elasticity(nk;
                                                                simplexify_model=true,
                                                                bcs_type=:mixed)
      Etag="elast"
      Ptag="poroelast"
      wallEDtag="wallED"
      wallPDtag="wallPD"
      wallENtag="wallEN"
      wallPNtag="wallPN"
      dirichlet_tags_Q=[wallPDtag,wallPNtag]
    elseif (problem==:vessel_bl || problem==:vessel_t)
      if (problem==:vessel_bl)
        model=generate_vessel_interstitium_model_block(:bl,nk,[:t,:r]) |> simplexify
      else
        model=generate_vessel_interstitium_model_block(:t,nk,[:b]) |> simplexify
      end
      Etag="vessel"
      Ptag="interstitium"
      wallEDtag="clamped_vessel"
      wallPDtag="clamped_interstitium"
      wallENtag="free_vessel"
      wallPNtag="free_interstitium"
      dirichlet_tags_Q=[wallPDtag,wallPNtag]
    end
  
    β_E, β_P, ξ_E, ξ_P, δ, θ1, θ2, Κ  = 
       scale_parameters(μ_E, μ_P, λ_E, λ_P, η, c_0, κ, α, 1.0; 
                                                             scale=scale, symmetric=true)

    # We may have here one more than one output
    # (as many as we like)
    output=Dict()
    rtol=1.0e-06
    niters,error_u,error_p,error_φ,ndofs=
       solve_stationary_hdiv_biot_elasticity(model,
                    order, β_E, β_P, ξ_E, ξ_P, θ1, θ2, δ, Κ, 1.0;
                    verbose=true,rtol=rtol,
                    β_U=β_U,
                    fe=fe,
                    solver=solver,
                    Etag=Etag,
                    Ptag=Ptag,
                    wallEDtag=wallEDtag,
                    wallPDtag=wallPDtag,
                    wallENtag=wallENtag,
                    wallPNtag=wallPNtag,
                    dirichlet_tags_Q=dirichlet_tags_Q)
    output["hh"]              = sqrt(2)/2^nk
    output["error_u"]         = error_u
    output["error_p"]         = error_p
    output["error_φ"]         = error_φ
    output["ndofs"]           = ndofs
    output["niters"]          = niters
    output["rtol"]            = rtol
    output
  end
end # module
