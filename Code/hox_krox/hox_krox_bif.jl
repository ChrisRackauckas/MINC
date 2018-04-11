using DiffEqBase, OrdinaryDiffEq
using Plots; plotly()
############### Parameters

#=
cₕ = 15.0
cₖ = 10.0
κₕ = .2
κₖ = 7.0
dₕ = 1.0
dₖ = 1.0
nₕ = 2
nₖ = 2
m = 2
=#

cₕ = 3
cₖ = 5
κₕ = 1
κₖ = 1
dₕ = 0.4
dₖ = 0.4
nₕ = 2
nₖ = 2
m = 2
aₕ = 0.05
aₖ = 0.05

#=
const cₕ = 3
const cₖ = 5
const κₕ = 1
const κₖ = 1
const dₕ = 0.4
const dₖ = 0.4
const nₕ = 2
const nₖ = 2
const m = 2
const aₕ = 0.05
const aₖ = 0.05
=#

################# Setup the problem

type HoxKrox
  RAin::Float64
end
function (p::HoxKrox)(du,u,p,t)
  gₕ,gₖ = u
  du[1] .= (cₕ.*gₕ.^nₕ .+ (κₕ*p.RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₕ*p.RAin).^m) - dₕ.*gₕ
  du[2] .= (cₖ.*gₖ.^nₖ .+ (κₖ*p.RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₖ*p.RAin).^m) - dₖ.*gₖ
end
u0 = [0.5,2.0]
hk = HoxKrox(2.0)
tspan = (0.0,500.0)
prob = ODEProblem(hk,u0,tspan)

####################### Setup Bif values

SShox = Vector{Float64}()
SSkrox = Vector{Float64}()
SShox2 = Vector{Float64}()
SSkrox2 = Vector{Float64}()
####################### Get steady states at different RAin

RAins = 0:0.01:2.5

for RAin_scalar in RAins
  prob.f.RAin = RAin_scalar
  sol = solve(prob,Tsit5(),save_everystep=false)
  push!(SShox,sol[end][1])
  push!(SSkrox,sol[end][2])
end

prob.u0 = [2.0,0.5]

for RAin_scalar in RAins
  prob.f.RAin = RAin_scalar
  sol = solve(prob,Tsit5(),save_everystep=false)
  push!(SShox2,sol[end][1])
  push!(SSkrox2,sol[end][2])
end

p1 = plot(RAins,[SShox,SSkrox],title="Steady State Krox Start",xaxis="RAin",
     yaxis="Steady State Value",label=["Hox" "Krox"],
     lw = 3);

p2 = plot(RAins,[SShox2,SSkrox2],title="Steady State Hox Start",xaxis="RAin",
     yaxis="Steady State Value",label=["Hox" "Krox"],
     lw = 3);
display(plot(p1,p2))

#=
####################### Get steady states at different RAin

RAin = ra_sol[end][:,:,2][:,21:61]
x = X[5,21:61]
RAins = RAin[5,:]
SShox = Vector{Float64}()
SSkrox = Vector{Float64}()
SShox2 = Vector{Float64}()
SSkrox2 = Vector{Float64}()

prob.u0 = [0.5,2.0]

for RAin_scalar in RAins
  prob.f.RAin = RAin_scalar
  sol = solve(prob,Tsit5(),save_everystep=false)
  push!(SShox,sol[end][1])
  push!(SSkrox,sol[end][2])
end

p21 = plot(x,[SShox,SSkrox],title="Steady State Hox/Krox vs RAin",xaxis="x (Anterior -> Posterior)",
     yaxis="Steady State Value",label=["Hox" "Krox"],
     lw = 3)

prob.u0 = [2.0,0.5]

for RAin_scalar in RAins
  prob.f.RAin = RAin_scalar
  sol = solve(prob,Tsit5(),save_everystep=false)
  push!(SShox2,sol[end][1])
  push!(SSkrox2,sol[end][2])
end

p22= plot(x,[SShox2,SSkrox2],title="Steady State Hox/Krox vs RAin",xaxis="x (Anterior -> Posterior)",
     yaxis="Steady State Value",label=["Hox" "Krox"],
     lw = 3)

display(plot(p21,p22))
=#
