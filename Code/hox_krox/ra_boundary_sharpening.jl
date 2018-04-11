using SpecialMatrices
using EllipsisNotation
using DiffEqBase, OrdinaryDiffEq, StochasticDiffEq
using Plots; pyplot()

############ Spatial Parameters

const x_start = -100
const x_end   = 400
const y_start = 0
const y_end   = 100
const dx = 5
const dy = 5
const N = Int((x_end-x_start) / dx + 1)
const M = Int((y_end-y_start) / dy +1)
const X = repmat(collect(x_start:dx:x_end)',M,1)
const Y = repmat(collect(y_start:dy:y_end),1,N)

############ RA Parameters

const DRA = 25.46
const β = 1
const kA = 0.002
const kdeg = 500 #1e5
#const kdeg = 500
const kmax = 1 #200
const γ = 2.5e-5
const f₀ = 400
const λ = 0.019
const n = 2
const xf = 400


# High noise

#const ϵout = 0.035
#const ϵin = 0.035

# Low noise

#const ϵout = 0.005
#const ϵin = 0.005

# No noise

const ϵout = 0.000
const ϵin = 0.000

############ Setup Diffusion Matrices

Ax = -full(Strang(N))
Ay = -full(Strang(M))
Ax[2,1] = 2
Ax[end-1,end] = 2
Ax[end,end] = -2*(1+dx*kA)
Ay[1,2] = 2
Ay[end,end-1] = 2
Ax ./= dx^2
Ay ./= dy^2

############# Spatial Constants

const Vmax = 1.0
V(x) = ifelse(x>xf-40,Vmax,zero(typeof(Vmax)))
const VRA = V.(X)
no_cyp(x) = ifelse(x<0 || (xf-40<x && x<=xf),true,false)
const NCYP_kmax = no_cyp.(X)

############# Hox-Krox Parameters

const cₕ = 7.5
const cₖ = 3.0
const κₕ = .4
const κₖ = 4.0
const dₕ = 0.4
const dₖ = 0.4
const nₕ = 2
const nₖ = 2
const m = 2
const aₕ = 0.15
const aₖ = 0.15
#const aₕ = 0.15
#const aₖ = 0.15

const tmp = zeros(M,N)
const Cyp = zeros(M,N)

u0 = zeros(M,N,2)
function ra_gradient(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]

  Cyp .= !NCYP_kmax.*(kdeg.*γ.*RAin.^n)./(1.+γ.*RAin.^n.+f₀.*exp.(-λ.*(xf.-X))) .+
         NCYP_kmax.*kmax
  dRAout .= DRA.*(RAout*Ax + Ay*RAout) .+ VRA .- (1.+β).*kA.*RAout .+ kA.*RAin
  dRAin  .= kA.*RAout .- (kA .+ Cyp).*RAin
  nothing
end

################ Solve ODE

tspan = (0.0,500.0)
prob = ODEProblem(ra_gradient,u0,tspan)
ra_sol = solve(prob,Tsit5(),save_everystep=false,force_dtmin=true)


RAin = ra_sol[end][:,:,2][:,21:61]
x = X[5,21:61]

################ Add Stochasticity

function ra_noise(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]

  dRAout .= ϵout .* RAout
  dRAin  .= ϵin  .* RAin
  nothing
end

u0 = ra_sol[end]
tspan2 = (0.0,50.0)
prob = SDEProblem(ra_gradient,ra_noise,u0,tspan2)
noisy_ra_sol = solve(prob,SRIW1(),save_everystep=false)

################ Plot Results

#const RAout = noisy_ra_sol[end][:,:,1]
#const RAin = noisy_ra_sol[end][:,:,2]
#surface(X,Y,RAout)
#surface(X,Y,RAin)
#surface(X,Y,Cyp)
#surface(X[:,21:61],Y[:,21:61],RAin[:,21:61],title="RAin Steady Distribution Gradient")
#contour(RAout,fill=true)
#RAout = sol[end][5,:,1]
#plot(X[1,:],RAout)

################# Hox-Krox Constant Initialization

u0 = 0.2ones(M,N,4)
u0[:,:,1] .= noisy_ra_sol[end][:,:,1]
u0[:,:,2] .= noisy_ra_sol[end][:,:,2]
u0[:,:,4] .*= 0.0

function hox_krox(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  gₕ = @view u[:,:,3]
  gₖ = @view u[:,:,4]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dgₕ = @view du[:,:,3]
  dgₖ = @view du[:,:,4]
  #Cyp .= !NCYP_kmax.*(kdeg.*γ.*RAin.^n)./(1.+γ.*RAin.^n.+f₀.*exp.(-λ.*(xf.-X))) .+
  #       NCYP_kmax.*kmax
  #dRAout .= DRA.*(RAout*Ax + Ay*RAout) .+ VRA .- (1.+β).*kA.*RAout .+ kA.*RAin
  #dRAin  .= kA.*RAout .- (kA .+ Cyp).*RAin
  dgₕ .= (cₕ.*gₕ.^nₕ .+ (κₕ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₕ*RAin).^m) - dₕ.*gₕ
  dgₖ .= (cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) - dₖ.*gₖ
  nothing
end

function hox_krox_noise(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  gₕ = @view u[:,:,3]
  gₖ = @view u[:,:,4]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dgₕ = @view du[:,:,3]
  dgₖ = @view du[:,:,4]

  #dRAout .= ϵout .* RAout
  #dRAin  .= ϵin  .* RAin
  dgₕ .= aₕ .* gₕ
  dgₖ .= aₖ .* gₖ
  nothing
end

tspan = (0.0,5000.0)
prob = ODEProblem(hox_krox,u0,tspan)
#prob = SDEProblem(hox_krox,hox_krox_noise,u0,tspan)
#sol = solve(prob,SRIW1(),save_everystep=false)
sol = solve(prob,Tsit5(),save_everystep=false)

gₕ = sol[end][:,:,3]
gₖ = sol[end][:,:,4]

#p1 = surface(X,Y,gₕ,title="hox")
#p2 = surface(X,Y,gₖ,title="krox")
#p1 = surface(X[:,21:61],Y[:,21:61],gₕ[:,21:61],title="hox")
#p2 = surface(X[:,21:61],Y[:,21:61],gₖ[:,21:61],title="krox")

#display(heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],title="Hoxb1a 2D Gradient, Constant Initialization",xlabel="Anterior -> Posterior",yaxis="Need to find out"))
#heatmap(X[1,21:61],Y[:,1],gₖ[:,21:61],title="Krox20 2D Gradient",
#        xlabel="Anterior -> Posterior",yaxis="Need to find out")

################# Hox-Krox

u0 = 0.2ones(M,N,4) + 0.05rand(M,N,4)
u0[:,:,1] .= noisy_ra_sol[end][:,:,1]
u0[:,:,2] .= noisy_ra_sol[end][:,:,2]
u0[:,:,4] .*= 0.0

function hox_krox(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  gₕ = @view u[:,:,3]
  gₖ = @view u[:,:,4]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dgₕ = @view du[:,:,3]
  dgₖ = @view du[:,:,4]
  #Cyp .= !NCYP_kmax.*(kdeg.*γ.*RAin.^n)./(1.+γ.*RAin.^n.+f₀.*exp.(-λ.*(xf.-X))) .+
  #       NCYP_kmax.*kmax
  #dRAout .= DRA.*(RAout*Ax + Ay*RAout) .+ VRA .- (1.+β).*kA.*RAout .+ kA.*RAin
  #dRAin  .= kA.*RAout .- (kA .+ Cyp).*RAin
  dgₕ .= (cₕ.*gₕ.^nₕ .+ (κₕ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₕ*RAin).^m) - dₕ.*gₕ
  dgₖ .= (cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) - dₖ.*gₖ
  nothing
end

function hox_krox_noise(du,u,p,t)
  RAout = @view u[:,:,1]
  RAin = @view u[:,:,2]
  gₕ = @view u[:,:,3]
  gₖ = @view u[:,:,4]
  dRAout = @view du[:,:,1]
  dRAin = @view du[:,:,2]
  dgₕ = @view du[:,:,3]
  dgₖ = @view du[:,:,4]

  #dRAout .= ϵout .* RAout
  #dRAin  .= ϵin  .* RAin
  dgₕ .= aₕ .* gₕ
  dgₖ .= aₖ .* gₖ
  nothing
end

################# Solve

tspan = (0.0,5000.0)
prob = ODEProblem(hox_krox,u0,tspan)
#prob = SDEProblem(hox_krox,hox_krox_noise,u0,tspan)
#sol = solve(prob,SRIW1(),save_everystep=false)
sol = solve(prob,Tsit5(),save_everystep=false)

gₕ = sol[end][:,:,3]
gₖ = sol[end][:,:,4]

#p1 = surface(X,Y,gₕ,title="hox")
#p2 = surface(X,Y,gₖ,title="krox")
#p1 = surface(X[:,21:61],Y[:,21:61],gₕ[:,21:61],title="hox")
#p2 = surface(X[:,21:61],Y[:,21:61],gₖ[:,21:61],title="krox")

p1 = heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],
        title="Hoxb1a 2D Gradient Before Sharpening",
        xlabel="Anterior/Posterior")
#heatmap(X[1,21:61],Y[:,1],gₖ[:,21:61],title="Krox20 2D Gradient",
#        xlabel="Anterior -> Posterior",yaxis="Need to find out")

########################### Hox Krox Noise

prob = SDEProblem(hox_krox,hox_krox_noise,sol[end],(0.0,40000.0))
sol2 = solve(prob,SRIW1(),save_everystep=false)

gₕ = sol2[end][:,:,3]
gₖ = sol2[end][:,:,4]

p2 = heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],
        title="Hoxb1a 2D Gradient After Sharpening",
        xlabel="Anterior/Posterior")
#heatmap(X[1,21:61],Y[:,1],gₖ[:,21:61],title="Krox20 2D Gradient",
#        xlabel="Anterior -> Posterior",yaxis="Need to find out")

display(plot(p1,p2,layout=grid(2,1)))
