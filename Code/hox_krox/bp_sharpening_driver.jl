#cd("D:\\OneDrive\\Current\\meanVar\\Code\\hox_krox")
cd("/home/crackauc/.julia/RAMeanVar/Code/hox_krox")
include("bp_sharpening_experiment_function.jl")
cd("data")
cd("noise_noise_plot")
cd("cyp_knocks")
############################# Solo Run

# 100.0 is sharp, 1.0 is unsharp
Vbp = 100.0#1.0
kdeg = 3.3e4
ϵₕₖ=0.30
sol = boundary_sharpening(Vbp=Vbp,ϵₕₖ=ϵₕₖ)
using JLD
#save("wildtype.jld", "u", sol.u, "t", sol.t)
save("bp_$(Int(Vbp))_epshk_$(ϵₕₖ).jld", "u", sol.u, "t", sol.t)

############################### Multi Run

Vbp = 5.0
ϵₕₖ = 0.3
using JLD
for Vbp in Vbps
  for i in 1:10
    sol = boundary_sharpening(Vbp=Vbp,ϵₕₖ=ϵₕₖ)
    save("bp_$(Vbp)_$i.jld", "u", sol.u, "t", sol.t)
  end
end

kdeg=3.3e4
using JLD
for i in 1:10
  sol = boundary_sharpening(kdeg=kdeg)
  save("cyp_$(kdeg)_$i.jld", "u", sol.u, "t", sol.t)
end


Vbps = 5.0:10.0:105.0
ϵₕₖ = 0.3

for Vbp in Vbps
    sol = boundary_sharpening(Vbp=Vbp,ϵₕₖ=ϵₕₖ)
    save("bp_$(Int(Vbp))_epshk_$(ϵₕₖ).jld", "u", sol.u, "t", sol.t)
end

Vbps = 5.0:10.0:105.0
ϵₕₖs = 0.1:0.025:0.3
for Vbp in Vbps, ϵₕₖ in ϵₕₖs
  sol = boundary_sharpening(Vbp=Vbp,ϵₕₖ=ϵₕₖ)
  save("bp_$(Int(Vbp))_epshk_$(ϵₕₖ).jld", "u", sol.u, "t", sol.t)
end
