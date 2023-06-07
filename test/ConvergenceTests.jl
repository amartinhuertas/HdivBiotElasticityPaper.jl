
# Julia script that can be used to reproduce the results reported in Section 7.1 of the paper 
# "Efficient and reliable divergence-conforming methods for an elasticity-poroelasticity interface problem" 
# by S.Badia, M. Hornkjöl, A. Khan, K.-A. Mardal, A. F. Martín, and R. Ruiz-Baier
#
using HdivBiotElasticityPaper
using Gridap 
using Printf

u_ex(x)  = VectorValue(sin(π*(x[1]+x[2])), cos(π*(x[1]^2+x[2]^2)))
p_ex(x)  = sin(π*x[1]+x[2])*sin(π*x[2])
φE_ex(x) = sin(2π*x[1]+x[2])*sin(4π*x[2])
φP_ex(x) = sin(8π*x[1]+x[2])*sin(10π*x[2])

function convergence_test(fe,nkmax;
                          u_ex=u_ex,
                          p_ex=p_ex,
                          φE_ex=φE_ex,
                          φP_ex=φP_ex,
                          simplexify_model=true,
                          parameter_combination=1,
                          β_U=(2.5*10^(2.0*(order-1)+1.0)),
                          order=1,
                          bcs_type=:full_dirichlet,
                          solver=:lu,
                          rtol=1.0e-06,
                          itmax=500)

    if (parameter_combination==1)
      κ=0.002
      α=1.0
      η=0.01
      μ_E = 0.37593984962406013
      μ_P = 33.557046979865774
      λ_E = 0.7297655904467051
      λ_P = 1644.2953020134212
      c_0 = 0.001
    elseif (parameter_combination==2)
      λ_P = 1.
      λ_E = 1.
      μ_P = 1.
      μ_E = 1.
      c_0 = 1.
      α = 1.
      κ = 1.
      η = 1.
    elseif (parameter_combination==3)
      # Exactly the same parameter-value combination as in Section 7.1 of manuscript
      λ_P = 2e+04
      λ_E = 1e+04
      μ_P = 10.
      μ_E = 20.
      c_0 = 1.
      α = 1.
      κ = 1.
      η = 1.
    elseif  (parameter_combination==4)
      # T-vessel application problem parameter values
      μ_E = 1.0e6
      μ_P = 1.0e3
      λ_E = 1.0e9
      λ_P = 1.0e6
      η   = 1.0
      c_0 = 1.0
      κ   = 1.0e-10
      α   = 1.0
    end

    β_E, β_P, ξ_E, ξ_P, δ, θ1, θ2, Κ  = 
                scale_parameters(μ_E, μ_P, λ_E, λ_P, η, c_0, κ, α, 1.0; 
                                                                    scale=true, 
                                                                    symmetric=true)

    eu   = Float64[]
    ep   = Float64[]
    ephi = Float64[]
    ru   = Float64[]
    rp   = Float64[]
    rphi = Float64[]
    hh   = Float64[]
    nn   = Int[]

    push!(ru,0.)
    push!(rp,0.)
    push!(rphi,0.)

    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model_unit_square_biot_elasticity(nk;
                                                         simplexify_model=simplexify_model,
                                                         bcs_type=bcs_type)
        push!(hh,sqrt(2)/2^nk)
        _,error_u,error_p,error_φ,ndofs=
              solve_stationary_hdiv_biot_elasticity(model,
                                      order,
                                      β_E, β_P, ξ_E, ξ_P, θ1, θ2, δ, Κ, 1.0;
                                      rtol=rtol,
                                      fe=fe,
                                      u_ex=u_ex,
                                      p_ex=p_ex,
                                      φE_ex=φE_ex,
                                      φP_ex=φP_ex,
                                      β_U=β_U,
                                      solver=solver,
                                      itmax=itmax)
        push!(nn,ndofs)
        println("******** Total DoFs: ", nn[nk])
        push!(eu,error_u)
        push!(ep,error_p)
        push!(ephi,error_φ)
        if nk>1
            push!(ru, log(eu[nk]/eu[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rp, log(ep[nk]/ep[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rphi, log(ephi[nk]/ephi[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("================================================================================")
    println("   DoF  &    h   &   e(u)   &  r(u)  &  e(p)  &  r(p)  &  e(phi)  &  r(phi)     ")
    println("================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f \n", nn[nk], hh[nk], eu[nk], ru[nk], ep[nk], rp[nk], ephi[nk], rphi[nk]);
    end
    println("================================================================================")
end

const parameter_combination=1  # 1, 2, 3, 4
const fe = :bdm                # :bdm, :raviart_thomas 
const nkmax = 5                # number of max uniform refinement levels 
const solver = :lu             # :lu, :pminres (Riesz mapping preconditioned)
const β_U=10.0                 # penalty parameter; Symmetric Interior Penalty DG
const simplexify_model=true    # true (triangles), false (quadrilaterals)
const order=1

# First touch Julia JIT compiler time takes a while, please be patient ...
# (Around 9 mins. with Julia 1.8.3 on a 12th Gen Intel(R) Core(TM) i7-1265U)
@time convergence_test(fe,nkmax;
                 simplexify_model=simplexify_model,
                 parameter_combination=parameter_combination,
                 β_U=β_U,
                 bcs_type=:mixed,
                 order=order,
                 solver=solver)