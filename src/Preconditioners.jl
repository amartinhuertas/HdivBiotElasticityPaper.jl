struct GridapLinearSolverPreconditioner{A,B}
  numerical_setup::A
  matrix::B
  function GridapLinearSolverPreconditioner(matrix)
    ls=LUSolver()
    ss=symbolic_setup(ls,matrix)
    ns=numerical_setup(ss,matrix)
    A=typeof(ns)
    B=typeof(matrix)
    new{A,B}(ns,matrix)
 end
end


function LinearAlgebra.mul!(z::AbstractVector{T},
              M::GridapLinearSolverPreconditioner,
              r::AbstractVector{T}) where T
  # Solve Mz=r ...
  Gridap.solve!(z,M.numerical_setup,r)
end
