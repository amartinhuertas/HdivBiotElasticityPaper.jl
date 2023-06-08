# Julia script that can be used to reproduce the results reported in Section 7.4 of the paper 
# "Efficient and reliable divergence-conforming methods for an elasticity-poroelasticity interface problem" 
# by S.Badia, M. Hornkjöl, A. Khan, K.-A. Mardal, A. F. Martín, and R. Ruiz-Baier
#

using HdivBiotElasticityPaper
using Gridap
using GridapGmsh

function main(;fem_formulation=:non_conforming,
               order=1,
               β_U=10.0,
               solution=:application,
               solver=:lu,
               u_ex=x->VectorValue(x[2],x[1]),
               p_ex=x->1.0,
               output_dir="output/",
               mesh_path="meshes/2D_T_vessel.msh",
               μ_E = 1.0e6,
               μ_P = 1.0e3,
               λ_E = 1.0e9,
               λ_P = 1.0e6,
               η   = 1.0,
               c_0 = 1.0,
               κ   = 1.0e-10,
               α   = 1.0,
               Δt = 0.01,
               T = 0.5,
               fe=:raviart_thomas,
               verbose=true)

  model = GmshDiscreteModel(mesh_path)
  set_boundary_vertices_entity_to_facets_entity!(model)

  k = fem_formulation==:non_conforming ? order : order+1
  θ = 0.5
  t = 0.0
  
  u_0 = VectorValue(0.0,0.0)
  p_0 = 0.0
  φ_0 = 0.0

  solve_transient_hdiv_biot_elasticity(
           model,
           k,t,T,Δt,θ,u_0,p_0,φ_0,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α;
           verbose=verbose,
           rtol=1.0e-06,
           fe=fe,
           β_U=β_U,
           solver=solver,
           fem_formulation=fem_formulation,
           solution=solution,
           u_ex=u_ex,
           p_ex=p_ex,
           output_dir=output_dir)
end


const fem_formulation=:non_conforming      # :non_conforming (BDM(k)-P(k)-P(k-1) or RT(k)-P(k)-P(k)), 
                                           # :conforming (P(k)-P(k-2)-P(k-1))
const fe = :bdm                            # :bdm, :raviart_thomas (only applies if fem_formulation=:non_conforming)
const β_U=10.0                             # penalty parameter; Symmetric Interior Penalty DG 
                                           # (only applies if fem_formulation=:non_conforming)
const solution=:application                # :application, :analytical
const solver = :lu                         # :lu, :pminres (Riesz mapping preconditioned)
const mesh_file = "meshes/2D_T_vessel.msh" # look at the "meshes/" dir for other mesh resolutions
const output_dir="output/"                 # directory of the output VTK files for later visualization in Paraview

# Physical parameters values
const μ_E = 1.0e6
const μ_P = 1.0e3
const λ_E = 1.0e9
const λ_P = 1.0e6
const η   = 1.0
const c_0 = 1.0
const κ   = 1.0e-10
const α   = 1.0

# Time parameters 
Δt  = 0.01  # time-step length 
T   = 0.5   # total simulation time 

# First touch Julia JIT compiler time takes a while, please be patient ...
# (Around 10 mins. with Julia 1.8.3 on a 12th Gen Intel(R) Core(TM) i7-1265U)
@time main(;
    fem_formulation=fem_formulation,
    fe=fe,
    β_U = β_U,
    solution=solution,
    solver=solver,
    mesh_path=mesh_file,
    output_dir=output_dir,
    μ_E = μ_E,
    μ_P = μ_P,
    λ_E = λ_E,
    λ_P = λ_P,
    η   = η,
    c_0 = c_0 ,
    κ   = κ,
    α   = α,
    Δt  = Δt,
    T   = T)
