"""
    struct PMinResSolver <: LinearSolver end

    Wrapper of the pminres solver available at Krylov.jl
"""
struct PMinResSolver{TP} <: LinearSolver
  P       :: TP # preconditioner
  itmax   :: Int
  verbose :: Int 
  rtol    :: Float64
  function PMinResSolver(P::TP;
                         itmax=500,
                         verbose=0,
                         rtol=1.0e-06) where TP 
    new{TP}(P,
            itmax,
            verbose,
            rtol)
  end                 
end 

struct PMinResSymbolicSetup{TP,TA} <: SymbolicSetup
  s::PMinResSolver{TP}
  A::TA # system matrix
end 

struct PMinResNumericalSetup{TP,TA} <: NumericalSetup
  s::PMinResSolver{TP} 
  A::TA # system matrix
end

Gridap.Algebra.symbolic_setup(solver::PMinResSolver, mat::AbstractMatrix) = PMinResSymbolicSetup(solver, mat)
Gridap.Algebra.numerical_setup(ss::PMinResSymbolicSetup, mat::AbstractMatrix) = PMinResNumericalSetup(ss.s, mat)

function Gridap.Algebra.numerical_setup!(ns::PMinResNumericalSetup, mat::AbstractMatrix)
  # Update preconditioner?
  #Gridap.Helpers.@notimplemented
  ns 
end

function Gridap.Algebra.numerical_setup!(ns::PMinResNumericalSetup, mat::SparseMatrixCSC)
  # Update preconditioner?
  #Gridap.Helpers.@notimplemented
  ns 
end

function custom_stopping_condition(solver::KrylovSolver, A, b, r, tol)
  mul!(r, A, solver.x)
  r .-= b                       # r := b - Ax
  bool = norm(r) ≤ tol*norm(b)  # tolerance based on the 2-norm of the residual
  @printf("||b-Ax||: %16.7e\n",norm(r))
  return bool
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::PMinResNumericalSetup,b::AbstractVector)
  r=copy(b)
  # println("norm(ns.A-ns.A')=", norm(ns.A-ns.A'))
  # println("max(ns.A-ns.A')=", maximum(ns.A-ns.A'))
  # println("norm(ns.A-ns.A')/norm(ns.A)=", norm(ns.A-ns.A')/norm(ns.A))
  minres_callback(solver) = custom_stopping_condition(solver, 
                                                      ns.A, 
                                                      b, 
                                                      r, 
                                                      ns.s.rtol)
  (xaux,stats)=minres(ns.A, b, M=ns.s.P; 
          itmax=ns.s.itmax, 
          verbose=ns.s.verbose,
          callback=minres_callback,
          atol=0.0,
          rtol=0.0,
          ratol=0.0,
          rrtol=0.0,
          etol=0.0)
  x .= xaux
end

"""
    struct PGmresSolver <: LinearSolver end

    Wrapper of the PGmres solver available at Krylov.jl
"""
struct PGmresSolver{TP} <: LinearSolver
  P       :: TP # preconditioner
  itmax   :: Int
  verbose :: Int 
  rtol    :: Float64
  function PGmresSolver(P::TP;
                         itmax=500,
                         verbose=0,
                         rtol=1.0e-06) where TP 
    new{TP}(P,
            itmax,
            verbose,
            rtol)
  end                 
end 

struct PGmresSymbolicSetup{TP,TA} <: SymbolicSetup
  s::PGmresSolver{TP}
  A::TA # system matrix
end 

struct PGmresNumericalSetup{TP,TA} <: NumericalSetup
  s::PGmresSolver{TP} 
  A::TA # system matrix
end

Gridap.Algebra.symbolic_setup(solver::PGmresSolver, mat::AbstractMatrix) = PGmresSymbolicSetup(solver, mat)
Gridap.Algebra.numerical_setup(ss::PGmresSymbolicSetup, mat::AbstractMatrix) = PGmresNumericalSetup(ss.s, mat)

function Gridap.Algebra.numerical_setup!(ns::PGmresNumericalSetup, mat::AbstractMatrix)
  # Update preconditioner?
  #Gridap.Helpers.@notimplemented
  ns 
end

function Gridap.Algebra.numerical_setup!(ns::PGmresNumericalSetup, mat::SparseMatrixCSC)
  # Update preconditioner?
  #Gridap.Helpers.@notimplemented
  ns 
end

function Gridap.Algebra.solve!(
  x::AbstractVector,ns::PGmresNumericalSetup,b::AbstractVector)
  r=copy(b)
  # println("norm(ns.A-ns.A')=", norm(ns.A-ns.A'))
  # println("max(ns.A-ns.A')=", maximum(ns.A-ns.A'))
  # println("norm(ns.A-ns.A')/norm(ns.A)=", norm(ns.A-ns.A')/norm(ns.A))
  gmres_callback(solver) = custom_stopping_condition(solver, 
                                                      ns.A, 
                                                      b, 
                                                      r, 
                                                      ns.s.rtol)
  (xaux,stats)=gmres(ns.A, b, M=ns.s.P; 
          itmax=ns.s.itmax, 
          verbose=ns.s.verbose,
          callback=gmres_callback,
          atol=0.0,
          rtol=ns.s.rtol)
  x .= xaux
end

