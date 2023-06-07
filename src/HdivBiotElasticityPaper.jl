module HdivBiotElasticityPaper
  using Gridap
  using Gridap.Algebra
  using SparseArrays
  using Printf
  using Krylov
  using LinearAlgebra
  include("DiscreteModels.jl")
  include("LinearSolvers.jl")
  include("Preconditioners.jl")  
  
  
  function compute_harmonic_μ_0(μ_E,μ_P)
    2.0*μ_E*μ_P/(μ_E+μ_P)
  end
  hmean(v,cf)      = ((v.⁺)*(cf.⁺/(cf.⁺ + cf.⁻)) + (v.⁻)*(cf.⁻/(cf.⁺ + cf.⁻)))
  hmean_swap(v,cf) = ((v.⁺)*(cf.⁻/(cf.⁺ + cf.⁻)) + (v.⁻)*(cf.⁺/(cf.⁺ + cf.⁻)))

  include("StationaryHdivBiotElasticityTools.jl")
  include("TransientHdivBiotElasticityTools.jl")
  
  export solve_stationary_hdiv_biot_elasticity
  export solve_transient_hdiv_biot_elasticity
  export generate_model_unit_square_biot_elasticity
  export scale_parameters
  export set_boundary_vertices_entity_to_facets_entity!
end
