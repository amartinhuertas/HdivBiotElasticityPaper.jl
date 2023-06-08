### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 38e1939a-3c91-4b7b-88fd-ae7e01eda602
# Run notebook in Pluto.jl "backward compatibility mode" 
# using the so-called "global environment" pattern
# see https://github.com/fonsp/Pluto.jl/wiki/%F0%9F%8E%81-Package-management 
# for more details
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()
end

# ╔═╡ 299c6e30-9107-4039-a11c-32ed6d9b1460
using DrWatson

# ╔═╡ 52a4a234-a99c-4e0d-b1e7-f54ddc11ff45
using DataFrames

# ╔═╡ 08a6ab1f-ef85-41a3-bb1a-933d578e4f93
using BSON

# ╔═╡ 696ffce0-e025-409c-8628-48b06fc7fd38
using PlutoUI

# ╔═╡ e82ad202-94eb-4f79-8aa3-7dd40375328b
using Plots

# ╔═╡ 2a36d1ee-af6b-424c-9ef7-8962c29fa910
using LaTeXStrings

# ╔═╡ 106ded38-479d-4961-a127-24c2b416aa9d
md"""## HdivBiotElasticityRieszPreconditioningTests
"""

# ╔═╡ 4eab9f96-8245-4d35-806f-21e34f75ae79
experiment_data_directory="HdivBiotElasticityRieszPreconditioningTests";

# ╔═╡ ebf65d95-b56a-4b69-ae99-127c3119d61b
commit_dirs=readdir(datadir(experiment_data_directory));

# ╔═╡ 7490a738-b0d1-4b9d-adb6-0fc27aadce2f
md"""
Select commit ID data directory: $(@bind commitID Select(commit_dirs)) 
"""

# ╔═╡ acbfd1ec-ef54-45fe-90e2-5bf9f919d728
df = collect_results(datadir(experiment_data_directory,commitID));

# ╔═╡ 602c6a8b-5a9d-4aba-9a73-c7af073feed0
cols_to_filter = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:ndofs,:niters,:solver,:scale];

# ╔═╡ 25889288-db3e-468d-8052-3ea6cfcb98f6
df_filtered = df[:,cols_to_filter];

# ╔═╡ 2f1b4812-7c6d-4510-b31c-234a7250e56e
df_filtered;

# ╔═╡ e4e03f80-42b4-459a-b75a-1e1364b3b659
df_filtered_cols = Dict(pairs(eachcol(df_filtered)));

# ╔═╡ 414e8521-662f-40f9-9c9e-def65d438fc7
params_possible_values=
	    Dict([k=>unique(df_filtered_cols[k]) 
			    for k in keys(df_filtered_cols)]);

# ╔═╡ 887ba2ff-c5b2-4fe3-ad20-03df4c9836b3
md"""
Select the parameter-value combination that you want to visualize!:

fe: $(@bind feval Select(params_possible_values[:fe])) 
βU: $(@bind βUval Select(params_possible_values[:β_U])) 
μE:   $(@bind μEval Select(params_possible_values[:μ_E]))
μP:   $(@bind μPval Select(params_possible_values[:μ_P]))
λE:   $(@bind λEval Select(params_possible_values[:λ_E]))
λP:   $(@bind λPval Select(params_possible_values[:λ_P]))
η:   $(@bind ηval Select(params_possible_values[:η]))
c0: $(@bind c0val Select(params_possible_values[:c_0]))
κ:   $(@bind κval Select(params_possible_values[:κ]))
α:   $(@bind αval Select(params_possible_values[:α])) 
solver: $(@bind solverval Select(params_possible_values[:solver];default=:pminres))
scale: $(@bind scaleval Select(params_possible_values[:scale];default=true))

Customize visualization

legend position: $(@bind lposition Select([:right, :left, :top, :bottom, :inside, :best, :legend, :topright, :topleft, :bottomleft, :bottomright];default=:topleft))

autoxlims: $(@bind autoxlims CheckBox(;default=true)) 
xlimleft: $(@bind xliml TextField((4,1);default="0.0"))
xlimright: $(@bind xlimr TextField((4,1);default="1.0"))

autoylims: $(@bind autoylims CheckBox(;default=false))
ylimbottom: $(@bind ylimb TextField((4,1);default="0.0"))
ylimtop: $(@bind ylimt TextField((4,1);default="500.0"))

logx: $(@bind logxval CheckBox(;default=true)) 
logy: $(@bind logyval CheckBox())

"""

# ╔═╡ 88d72926-df2f-4a47-9de9-39414f079498
ffilter=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,solver,scale)->(fe==feval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval && solver==solverval && scale==scaleval);

# ╔═╡ 5748a503-4534-462a-b485-386d68fe857d
function generate_labels(params_possible_values)
  count=Dict{Any,Int}()
  for (key,val) in zip(keys(params_possible_values),values(params_possible_values))
	  count[key]=length(val)
  end
  dl=dict_list(params_possible_values)
  labels=Vector{String}(undef,length(dl))
  for (i,d) in enumerate(dl)
	label=""
	for (key,val) in d
	  if key !== :fe && count[key]>1 
		if (key==:β_U)
		  label=label * L" \beta_{\mathbf{u}}" * "=$(val)"
		elseif (key==:nk)
		  label=label * L" \ell" * "=$(val)" 	
		else		
	      label=label * " $(key)=$(val)"
		end
	  end 
	end
	labels[i]=label  
  end
  labels
end 

# ╔═╡ 2c7a1e9b-ebe7-4631-a14d-0c2ea0d7293b
function get_x_y(xparam, yparam, ffilter, df)
  df_filtered = filter(ffilter,df)
  df_filtered_cols = Dict(pairs(eachcol(df_filtered)))
  @assert xparam in keys(df_filtered_cols)
  @assert yparam in keys(df_filtered_cols)
	
  params_possible_values=
	    Dict([k=>unique(df_filtered_cols[k]) 
			    for k in keys(df_filtered_cols) if k != xparam && k != yparam])
	
  dl=dict_list(params_possible_values)

  # The following code is general enough so that for 
  # fixed (xparam, yparam) there might be several 
  # possible combinations for the rest of parameters
  # after applying ffilter. In such a case we generate
  # as many curves as combinations of the rest of 
  # parameter values.
  xall=[]
  yall=[]	
  for d in dl
      function f(a...)
		  all(a .== values(d))
	  end 
	  ffilter_current_d = collect(keys(d))=>f
	  df_tmp=filter(ffilter_current_d,df_filtered)
	  sort!(df_tmp,[xparam,])
      x = df_tmp[!,xparam]
      y = df_tmp[!,yparam]
      push!(xall,x)
	  push!(yall,y)
  end
  (xall,yall,params_possible_values)
end


# ╔═╡ 9478d88e-bc39-40b2-be09-3789e4fff51e
function plot_xparam_versus_yparam(xparam,yparam,xaxis,yaxis,ffilter,df;
                                   autoxlims=true,
                                   autoylims=true,
                                   xliml=0.0,
                                   xlimr=1.0,
                                   ylimb=0.0,
                                   ylimt=500.0,
                                   xlabel="$xparam",
								   ylabel="$yparam")
  f = plot()
  x,y,params_possible_values = get_x_y(xparam,yparam,ffilter,df)
  labels=generate_labels(params_possible_values)
  @assert length(x)==length(y)
  @assert length(labels)==length(y)	
  for (xi,yi,li) in zip(x,y,labels)
    plot!(xi,yi,xaxis=xaxis,yaxis=yaxis,label=li,markershape=:auto)
  end 
  plot!(xlabel=xlabel,ylabel=ylabel,legend=lposition)
  if (!autoxlims)
    xlims!((xliml,xlimr))
  end	  
  if (!autoylims)
    ylims!((ylimb,ylimt))
  end	
  f
end

# ╔═╡ 949a60c0-6be3-4e0e-8f87-264276e0943a
begin
	iters_versus_dofs=plot_xparam_versus_yparam(:ndofs,:niters,
		                      (logxval ? :log10 : :none),
		                      (logyval ? :log10 : :none),
	                          ffilter,
		                      df_filtered;
	                          autoxlims=autoxlims,
							  autoylims=false,
	                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
							  ylimb=0,ylimt=500,
							  ylabel=L"\#\mathrm{iterations}",
	                          xlabel=L"\mathrm{DoF}",)
	savefig(iters_versus_dofs, "iters_vs_dofs.pdf");
	iters_versus_dofs
end

# ╔═╡ f899fd95-9069-4fc1-92a7-db567d580e25
begin 
  cols_to_filter_beta = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:niters,:solver,:scale];
  ffilter_beta=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,nk,solver,scale)->(fe==feval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  && solver==:pminres && scale==scaleval);
  df_filtered_beta = df[:,cols_to_filter_beta];
  iters_vs_beta=plot_xparam_versus_yparam(:β_U,:niters,
	                      (logxval ? :log10 : :none),
	                      (logyval ? :log10 : :none),
                          ffilter_beta,
	                      df_filtered_beta;
                          autoxlims=true,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=0,ylimt=500,						  
						  ylabel=L"\#\mathrm{iterations}",
	                      xlabel=L"\beta_{\mathbf{u}}",
)
savefig(iters_vs_beta, "iters_vs_beta.pdf");
	iters_vs_beta
end	

# ╔═╡ 4d27e8ec-c25e-4643-bc91-e7124ef92ccb
begin 
  cols_to_filter_beta_error = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:error_u,:solver,:scale];
  ffilter_beta_error=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,nk,solver,scale)->(fe==feval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  &&  scale==scaleval);
  df_filtered_beta_error = df[:,cols_to_filter_beta_error];
  error_u_vs_beta=plot_xparam_versus_yparam(:β_U,:error_u,
	                      (logxval ? :log10 : :none),
	                      :log10,
                          ffilter_beta_error,
	                      df_filtered_beta_error;
                          autoxlims=autoxlims,
						  autoylims=true,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-04,ylimt=1.0e+01,
            ylabel=L"\|\|\| \mathbf{u}-\mathbf{u}_h \|\|\|_{\ast \ast}",
            xlabel=L"\beta_{\mathbf{u}}")	
	savefig(error_u_vs_beta, "error_u_vs_beta.pdf");
	error_u_vs_beta
end

# ╔═╡ 472fddbf-591d-4275-ab2c-1bffc9f029fc
begin 
  cols_to_filter_beta_error_p = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:error_p,:solver,:scale];
  ffilter_beta_error_p=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,nk,solver,scale)->(fe==feval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  &&  scale==scaleval);
  df_filtered_beta_error_p = df[:,cols_to_filter_beta_error_p];
  error_p_vs_beta=plot_xparam_versus_yparam(:β_U,:error_p,
	                      (logxval ? :log10 : :none),
	                      :log10,
                          ffilter_beta_error_p,
	                      df_filtered_beta_error_p;
                          autoxlims=autoxlims,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-05,ylimt=1.0e-01,
            ylabel=L"\|\|\| p^{\mathrm{P}}-p^{\mathrm{P}}_h \|\|\|_{\ast \ast}",
            xlabel=L"\beta_{\mathbf{u}}")	
		savefig(error_p_vs_beta, "error_p_vs_beta.pdf");
	error_p_vs_beta
end

# ╔═╡ 1ce8eaab-112a-417b-9a0e-2bd0be97d2c1
begin 
  cols_to_filter_beta_error_phi = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:error_φ,:solver,:scale];
  ffilter_beta_error_phi=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:nk,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,nk,solver,scale)->(fe==feval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  &&  scale==scaleval);
  df_filtered_beta_error_phi = df[:,cols_to_filter_beta_error_phi];
  error_phi_vs_beta=plot_xparam_versus_yparam(:β_U,:error_φ,
	                      (logxval ? :log10 : :none),
	                      :log10,
                          ffilter_beta_error_phi,
	                      df_filtered_beta_error_phi;
                          autoxlims=autoxlims,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-02,ylimt=1.0e+03,
            ylabel=L"\|\|\| \varphi-\varphi_h \|\|\|_{\ast \ast}",
            xlabel=L"\beta_{\mathbf{u}}")
	savefig(error_phi_vs_beta, "error_phi_vs_beta.pdf");
	error_phi_vs_beta
end

# ╔═╡ 93686b31-35a2-41ff-874d-f4f76c212509
begin 
  cols_to_filter_error = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:ndofs,:error_u, :solver,:scale];
  ffilter_error=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,solver,scale)->(fe==feval && β_U==βUval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval && scale==scaleval);
  df_filtered_error = df[:,cols_to_filter_error]
  plot_xparam_versus_yparam(:ndofs,:error_u,
	                      :log10,
	                      :log10,
                          ffilter_error,
	                      df_filtered_error;
                          autoxlims=true,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-04,ylimt=1.0e-02,
        ylabel=L"\|\|\| \mathbf{u}-\mathbf{u}_h \|\|\|_{\ast \ast}",
	    xlabel=L"\mathrm{DoF}",
   )	
end

# ╔═╡ c1f5b97d-e391-45de-bf48-ffe7cb9e53a9
begin 
  cols_to_filter_error_p = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:ndofs,:error_p,:solver,:scale];
  ffilter_error_p=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,solver,scale)->(fe==feval && β_U==βUval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  && scale==scaleval);
  df_filtered_error_p = df[:,cols_to_filter_error_p]
plot_xparam_versus_yparam(:ndofs,:error_p,
	                      :log10,
	                      :log10,
                          ffilter_error_p,
	                      df_filtered_error_p;
                          autoxlims=true,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-05,ylimt=1.0e-02,
        ylabel=L"\|\|\| p^{\mathrm{P}}-p^{\mathrm{P}}_h \|\|\|_{\ast \ast}",
	    xlabel=L"\mathrm{DoF}",)	
end

# ╔═╡ eabba61d-cc2b-41cb-8e69-6b7f24591c26
begin 
  cols_to_filter_error_φ = [:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:ndofs,:error_φ,:solver,:scale];
  ffilter_error_φ=[:fe,:β_U,:μ_E,:μ_P,:λ_E,:λ_P,:η,:c_0,:κ,:α,:solver,:scale]=>(fe,β_U,μ_E,μ_P,λ_E,λ_P,η,c_0,κ,α,solver,scale)->(fe==feval && β_U==βUval && μ_E==μEval && μ_P==μPval && λ_E==λEval && λ_P==λPval && η==ηval && c_0==c0val && κ==κval && α==αval  && scale==scaleval);
  df_filtered_error_φ = df[:,cols_to_filter_error_φ]
plot_xparam_versus_yparam(:ndofs,:error_φ,
	                      :log10,
	                      :log10,
                          ffilter_error_φ,
	                      df_filtered_error_φ;
                          autoxlims=true,
						  autoylims=false,
                          xliml=parse(Float64, xliml),xlimr=parse(Float64, xlimr),
						  ylimb=1.0e-02,ylimt=1.0e+00,
						  ylabel=L"\|\|\| \varphi-\varphi_h \|\|\|_{\ast \ast}",
	    xlabel=L"\mathrm{DoF}",)	
end

# ╔═╡ Cell order:
# ╠═106ded38-479d-4961-a127-24c2b416aa9d
# ╠═38e1939a-3c91-4b7b-88fd-ae7e01eda602
# ╠═299c6e30-9107-4039-a11c-32ed6d9b1460
# ╠═52a4a234-a99c-4e0d-b1e7-f54ddc11ff45
# ╠═08a6ab1f-ef85-41a3-bb1a-933d578e4f93
# ╠═696ffce0-e025-409c-8628-48b06fc7fd38
# ╠═e82ad202-94eb-4f79-8aa3-7dd40375328b
# ╠═2a36d1ee-af6b-424c-9ef7-8962c29fa910
# ╠═4eab9f96-8245-4d35-806f-21e34f75ae79
# ╠═ebf65d95-b56a-4b69-ae99-127c3119d61b
# ╠═7490a738-b0d1-4b9d-adb6-0fc27aadce2f
# ╠═acbfd1ec-ef54-45fe-90e2-5bf9f919d728
# ╠═602c6a8b-5a9d-4aba-9a73-c7af073feed0
# ╠═25889288-db3e-468d-8052-3ea6cfcb98f6
# ╠═2f1b4812-7c6d-4510-b31c-234a7250e56e
# ╠═e4e03f80-42b4-459a-b75a-1e1364b3b659
# ╠═414e8521-662f-40f9-9c9e-def65d438fc7
# ╠═88d72926-df2f-4a47-9de9-39414f079498
# ╠═949a60c0-6be3-4e0e-8f87-264276e0943a
# ╠═f899fd95-9069-4fc1-92a7-db567d580e25
# ╠═4d27e8ec-c25e-4643-bc91-e7124ef92ccb
# ╠═472fddbf-591d-4275-ab2c-1bffc9f029fc
# ╠═1ce8eaab-112a-417b-9a0e-2bd0be97d2c1
# ╠═887ba2ff-c5b2-4fe3-ad20-03df4c9836b3
# ╠═93686b31-35a2-41ff-874d-f4f76c212509
# ╠═c1f5b97d-e391-45de-bf48-ffe7cb9e53a9
# ╠═eabba61d-cc2b-41cb-8e69-6b7f24591c26
# ╠═9478d88e-bc39-40b2-be09-3789e4fff51e
# ╠═5748a503-4534-462a-b485-386d68fe857d
# ╟─2c7a1e9b-ebe7-4631-a14d-0c2ea0d7293b
