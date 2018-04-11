using SpecialMatrices
using EllipsisNotation
using DiffEqBase, OrdinaryDiffEq, StochasticDiffEq
using Plots; plotly()

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
const β = .017 # control steepness, higher is steeper
const kA = 0.002
const kdeg = 1e5
const kmax = 1
const n = 2
const xf = 400

const kp = 0.1
const gamma = 100
const mon = 3.0
const moff = 0.0013
const rdeg1 = 0.0001
const rdeg2 = 0.0001
const jalpha = .85
const jbeta = .850
const bpdeg1 = 0.01
const d = 0.1
const e = 1
# 100.0 is sharp, 1.0 is unsharp
const Vbp = 1.0#1.0
const Vr = 2e-2

# Low noise

const ϵout = 0.03
const ϵin = 0.03
const ϵR = 0.03

############ Setup Diffusion Matrices

const Ax = -full(Strang(N))
const Ay = -full(Strang(M))
Ax[2,1] = 2
Ax[end-1,end] = 2
Ax[end,end] = -2*(1+dx*kA)
Ay[1,2] = 2
Ay[end,end-1] = 2
Ax ./= dx^2
Ay ./= dy^2

const Dxtmp = zeros(M,N)
const Dytmp = zeros(M,N)

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
# 10.75
const aₕ = .2
const aₖ = .2

u0 = zeros(M,N,6)

function ra_gradient(du,u,p,t) # no matrix cache for autodiff, better for stiff
  RAout = view(u,:,:,1)
  RAin = view(u,:,:,2)
  R = view(u,:,:,3)
  RAR = view(u,:,:,4)
  BP = view(u,:,:,5)
  RABP = view(u,:,:,6)
  dRAout = view(du,:,:,1)
  dRAin = view(du,:,:,2)
  dR = view(du,:,:,3)
  dRAR = view(du,:,:,4)
  dBP = view(du,:,:,5)
  dRABP = view(du,:,:,6)

  # Doesn't matter because all the time is in the factorization anyways
  Dxtmp = RAout*Ax
  Dytmp = Ay*RAout

  Threads.@threads for i in eachindex(dRAout)
    @inbounds diffRA     = DRA*(Dxtmp[i] + Dytmp[i])
    @inbounds dRAout[i]  = VRA[i] - β*RAout[i] + kp*RAin[i] + diffRA
    @inbounds dRAin[i]   = β*RAout[i] - kp*RAin[i] - kdeg*(RAR[i]/(gamma+RAR[i]))*RAin[i] - mon*RAin[i]*BP[i] + moff*RABP[i] - rdeg1*RAin[i]
    @inbounds dR[i]      = Vr - rdeg2*R[i]- jalpha*RABP[i]*R[i] + jbeta*BP[i]*RAR[i]
    @inbounds dRAR[i]    = jalpha*RABP[i]*R[i] - jbeta*BP[i]*RAR[i]
    @inbounds dBP[i]     = Vbp - bpdeg1*BP[i] - mon*RAin[i]*BP[i] + moff*RABP[i] + jalpha*RABP[i]*R[i] - jbeta*BP[i]*RAR[i] + ((d*RAR[i])/(e+RAR[i]))
    @inbounds dRABP[i]   = mon*RAin[i]*BP[i] - moff*RABP[i] - jalpha*RABP[i]*R[i] + jbeta*BP[i]*RAR[i]
  end
  nothing
end

function ra_gradient2(du,u,p,t) # better for non-stiff
  RAout = view(u,:,:,1)
  RAin = view(u,:,:,2)
  R = view(u,:,:,3)
  RAR = view(u,:,:,4)
  BP = view(u,:,:,5)
  RABP = view(u,:,:,6)
  dRAout = view(du,:,:,1)
  dRAin = view(du,:,:,2)
  dR = view(du,:,:,3)
  dRAR = view(du,:,:,4)
  dBP = view(du,:,:,5)
  dRABP = view(du,:,:,6)

  A_mul_B!(Dxtmp,RAout,Ax)
  A_mul_B!(Dytmp,Ay,RAout)

  Threads.@threads for i in eachindex(dRAout)
    @inbounds diffRA     = DRA*(Dxtmp[i] + Dytmp[i])
    @inbounds dRAout[i]  = VRA[i] - β*RAout[i] + kp*RAin[i] + diffRA
    @inbounds dRAin[i]   = β*RAout[i] - kp*RAin[i] - kdeg*(RAR[i]/(gamma+RAR[i]))*RAin[i] - mon*RAin[i]*BP[i] + moff*RABP[i] - rdeg1*RAin[i]
    @inbounds dR[i]      = Vr - rdeg2*R[i]- jalpha*RABP[i]*R[i] + jbeta*BP[i]*RAR[i]
    @inbounds dRAR[i]    = jalpha*RABP[i]*R[i] - jbeta*BP[i]*RAR[i]
    @inbounds dBP[i]     = Vbp - bpdeg1*BP[i] - mon*RAin[i]*BP[i] + moff*RABP[i] + jalpha*RABP[i]*R[i] - jbeta*BP[i]*RAR[i] + ((d*RAR[i])/(e+RAR[i]))
    @inbounds dRABP[i]   = mon*RAin[i]*BP[i] - moff*RABP[i] - jalpha*RABP[i]*R[i] + jbeta*BP[i]*RAR[i]
  end
  nothing
end

################ Solve ODE

println("Solve the ODE")
tspan = (0.0,1e9)
prob = ODEProblem(ra_gradient,u0,tspan)
@time ra_sol = solve(prob,Rosenbrock23(),save_everystep=false,dt=0.5,
                     progress=true,abstol=1e-3,reltol=1e-1,progress_steps=1)

################ Add Stochasticity

println("Solve the gradient SDE")
function ra_noise(du,u,p,t)
  RAout = view(u,:,:,1)
  RAin = view(u,:,:,2)
  RAR = view(u,:,:,4)
  dRAout = view(du,:,:,1)
  dRAin = view(du,:,:,2)
  dRAR = view(du,:,:,4)

  Threads.@threads for i in eachindex(dRAin)
    @inbounds dRAout[i] = ϵout * RAout[i]
    @inbounds dRAin[i] = ϵin * RAin[i]
    @inbounds dRAR[i] = ϵR * RAR[i]
  end
  nothing
end

#################

u0 = ra_sol[end]
tspan2 = (0.0,500.0)
prob = SDEProblem(ra_gradient2,ra_noise,u0,tspan2)
@time noisy_ra_sol = solve(prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),
                           save_everystep=false,progress_steps=10_000,
                           progress=true,abstol=1e-2,reltol=1e-2)

################ Plot Results

surface(X[:,21:61],Y[:,21:61],noisy_ra_sol[end][:,21:61,1],title="RAout")
surface(X[:,21:61],Y[:,21:61],noisy_ra_sol[end][:,21:61,2],title="RAin")
surface(X[:,21:61],Y[:,21:61],noisy_ra_sol[end][:,21:61,4],title="RAR")
surface(X[:,21:61],Y[:,21:61],noisy_ra_sol[end][:,21:61,5],title="BP")


################# Hox-Krox DiffEqs

function hox_krox(du,u,p,t)
  RAR = view(u,:,:,4)
  gₕ = view(u,:,:,7)
  gₖ = view(u,:,:,8)
  dgₕ = view(du,:,:,7)
  dgₖ = view(du,:,:,8)

  Threads.@threads for i in eachindex(dgₕ)
    @inbounds dgₕ[i] = (cₕ*gₕ[i]^nₕ + (κₕ*RAR[i])^m) / (1 + cₕ*gₕ[i]^nₕ + cₖ*gₖ[i]^nₖ + (κₕ*RAR[i])^m) - dₕ*gₕ[i]
    @inbounds dgₖ[i] = (cₖ*gₖ[i]^nₖ + (κₖ*RAR[i])^m) / (1 + cₕ*gₕ[i]^nₕ + cₖ*gₖ[i]^nₖ + (κₖ*RAR[i])^m) - dₖ*gₖ[i]
  end
  nothing
end

function hox_krox_noise(du,u,p,t)
  RAout = view(u,:,:,1)
  RAin = view(u,:,:,2)
  RAR = view(u,:,:,4)
  gₕ = view(u,:,:,7)
  gₖ = view(u,:,:,8)
  dRAout = view(du,:,:,1)
  dRAin = view(du,:,:,2)
  dRAR = view(du,:,:,4)
  dgₕ = view(du,:,:,7)
  dgₖ = view(du,:,:,8)

  Threads.@threads for i in eachindex(dRAin)
    @inbounds dgₕ[i] = aₕ[i] * gₕ[i]
    @inbounds dgₖ[i] = aₖ[i] * gₖ[i]
  end
  nothing
end

################# Hox-Krox

u0 = 0.16ones(M,N,8) + 0.20rand(M,N,8) + 0.0005
u0[:,:,1] .= noisy_ra_sol[end][:,:,1]
u0[:,:,2] .= noisy_ra_sol[end][:,:,2]
u0[:,:,3] .= noisy_ra_sol[end][:,:,3]
u0[:,:,4] .= noisy_ra_sol[end][:,:,4]
u0[:,:,5] .= noisy_ra_sol[end][:,:,5]
u0[:,:,6] .= noisy_ra_sol[end][:,:,6]
u0[:,1:46,7] = ones(X[:,1:46])
u0[:,:,8] .*= 0.0

println("Solve the Hox-Krox ODE")

tspan = (0.0,5000.0)
prob = ODEProblem(hox_krox,u0,tspan)

@time sol = solve(prob,Tsit5(),save_everystep=false,progress=true,
            abstol=1e-3,reltol=1e-2,progress_steps=100)

gₕ = sol[end][:,:,7]

p1 = heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],
        title="Hoxb1a 2D Gradient Before Sharpening",
        xlabel="Anterior/Posterior")

########################### Hox Krox Noise

println("Solve the Hox-Krox SDE")

prob = SDEProblem(hox_krox,hox_krox_noise,sol[end],(0.0,10000.0))
@time sol2 = solve(prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),
      progress_steps=1000,saveat=100.0,
      progress=true,abstol=1e-2,reltol=1e-2)

gₕ = sol2[end][:,:,7]

p2 = heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],
        title="Hoxb1a 2D Gradient After Sharpening",
        xlabel="Anterior/Posterior")

display(plot(p1,p2,layout=grid(2,1)))


### Intermediate Checks

#surface(X[:,21:61],Y[:,21:61],ra_sol[end][:,21:61,1],title="RAout")
#surface(X[:,21:61],Y[:,21:61],ra_sol[end][:,21:61,2],title="RAin")
#surface(X[:,21:61],Y[:,21:61],ra_sol[end][:,21:61,4],title="RAR")
#
#u = copy(ra_sol[end])
#du = similar(u)
#t = 0
#@time ra_gradient(du,u,p,t)
#maximum(du)

################# Hox-Krox Constant Initialization

#println("Solve the Hox-Krox ODE, Constant Initialization")
#noisy_ra_sol = ra_sol
#u0 = 0.2ones(M,N,8)
#u0[:,:,1] .= noisy_ra_sol[end][:,:,1]
#u0[:,:,2] .= noisy_ra_sol[end][:,:,2]
#u0[:,:,3] .= noisy_ra_sol[end][:,:,3]
#u0[:,:,4] .= noisy_ra_sol[end][:,:,4]
#u0[:,:,5] .= noisy_ra_sol[end][:,:,5]
#u0[:,:,6] .= noisy_ra_sol[end][:,:,6]
#u0[:,:,8] .*= 0.0
#tspan = (0.0,5000.0)
#prob = ODEProblem(hox_krox,u0,tspan)
#sol = solve(prob,Tsit5(),save_everystep=false,progress=true,abstol=1e-3,reltol=1e-1)
#gₕ = sol[end][:,:,7]
#display(heatmap(X[1,31:61],Y[:,1],gₕ[:,31:61],title="Hoxb1a 2D Gradient, Constant Initialization",xlabel="Anterior -> Posterior",yaxis="Need to find out"))
