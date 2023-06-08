
using DrWatson

# The following macro call let us execute the script with
# the project's environment even if we ran the julia REPL
# without the --project=... flag
@quickactivate "HdivBiotElasticityPaper"

function replace_strings_by_symbols(r)
  d = Dict{Symbol,Any}()
  for (k,v) in r
    d[Symbol(k)] = v
  end
  d
end

# Define parameter-value combinations for the experiment.
# Parameter-value combinations s.t. corresponding results
# are already available in the data folder are not re-run.
# You have to eliminate them to be re-run
function generate_param_dicts()
   params = Dict{Symbol,Any}(
     :problem => [:vessel_t],
     :solver  => [:lu,:pminres],
     :fe      => [:bdm],
     :nk      => [2,4,8],
     :order   => [1],
     :β_U     => [10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0,5000.0],
     :μ_E     => [1.0e6],
     :μ_P     => [1.0e3],
     :λ_E     => [1.0e9],
     :λ_P     => [1.0e6],
     :η       => [1.],
     :c_0     => [1.0],
     :κ       => [1.0e-10],
     :α       => [1.0],
     :scale   => [true,false],
   )
   dict_list(params)
end

# Defines the Driver module with the driver(...) function inside
# The computational heavy stuff and the actual code of the experiment
# at hand is here.
include("driver.jl")

function run_experiment(params)
  outfile = datadir("HdivBiotElasticityRieszPreconditioningTests",gitdescribe(),
                    savename("HdivBiotElasticityRieszPreconditioningTests",params,"bson"))
  if isfile(outfile)
    println("$outfile (done already)")
    return nothing
  end
  println("$outfile (running)")
  @unpack problem, solver,fe,nk,order,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,scale = params
  dict = Driver.driver(problem,solver,fe,nk,order,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,scale)
  merge!(dict,params)
  # @tagsave: add current git commit to dict and then save.
  # "replace_strings_by_symbols" is required to ensure that
  # the dictionary is not of type Dict{Any} but of type
  # Dict{Symbol}
  @tagsave(outfile,replace_strings_by_symbols(dict))
  println(" (done)")
end

# Run all parameter value combinations
dicts=generate_param_dicts()
for params in dicts
  GC.gc()
  run_experiment(params)
end
