cd("D:\\OneDrive\\Current\\meanVar\\Code\\hox_krox")
include("bp_sharpening_experiment_function.jl")
cd("D:\\OneDrive\\Backups\\ArchivedProjects\\UCI\\MeanVar Independance\\arxdata\\hox_krox_data")
using JLD

###################### Analysis Consts

const x = X[1,31:61]-50
const y = Y[:,1]
const boundary_xs = Vector{Float64}(y)

function calculate_boundary_stds(u)
  boundary_stds = Vector{Float64}(length(u))
  for ti in eachindex(u)
    gₕ = u[ti][:,31:61,7]
    for i in 1:size(gₕ,1)
      cur_g = gₕ[i,:] .< 0.5
      idx = findfirst(cur_g)
      boundary_xs[i] = x[idx]
    end
    boundary_stds[ti] = sqrt(var(boundary_xs))
  end
  boundary_stds
end

function calculate_boundary_stds2(u)
  boundary_stds = Vector{Float64}(length(u))
  for ti in eachindex(u)
    bad_distances = Vector{Float64}(0)
    gₕ = u[ti][:,31:61,7]
    numy = size(gₕ,1)
    numx = size(gₕ,2)
    for i in 1:numy
      cur_g = gₕ[i,:] .< 0.5
      idx = findfirst(cur_g)
      idx2 = findlast(.!(cur_g))
      mid_idx = (idx+idx2)/2
      midpoint = (x[idx] + x[idx2])/2
      bad_left_idxs  = (1:round(Int,mid_idx,RoundDown))[cur_g[1:round(Int,mid_idx,RoundDown)]]
      bad_right_idxs = (round(Int,mid_idx,RoundUp):numx)[.!cur_g[round(Int,mid_idx,RoundUp):numx]]
      bad_left_distances = x[bad_left_idxs] - midpoint
      bad_right_distances = x[bad_right_idxs] - midpoint
      append!(bad_distances,[bad_left_distances;bad_right_distances])
    end
    boundary_stds[ti] = 2std(bad_distances)/numy
  end
  boundary_stds
end

function calculate_boundary_stds3(u,reduction=mean)
  boundary_stds = Vector{Float64}(length(u))
  for ti in eachindex(u)
    gₕ = u[ti][:,31:61,7]
    for i in 1:size(gₕ,1)
      cur_g = gₕ[i,:] .< 0.5
      idx = findfirst(cur_g)
      idx2 = findlast(.!(cur_g))
      boundary_xs[i] = max(x[idx2] - x[idx],0)
    end
    boundary_stds[ti] = reduction(boundary_xs)
  end
  boundary_stds
end

function calculate_boundary_stds4(u,cell_length=10)
  boundary_stds = Vector{Float64}(length(u))
  for ti in eachindex(u)
    gₕ = u[ti][:,31:61,7]
    for i in 1:size(gₕ,1)
      cur_g = gₕ[i,:] .< 0.5
      idx = findfirst(cur_g)
      idx2 = findlast(.!(cur_g))
      boundary_xs[i] = max(x[idx2] - x[idx],0)
    end
    boundary_stds[ti] = sum(boundary_xs .> cell_length)
  end
  boundary_stds
end


effective_σRA(u) = var(u[end][:,46,4])

####################### Load the data

using JLD
data1 = load("wildtype/bp_100.0_10.jld")
data2 = load("bp_knocks/bp_5.0_7.jld")
data3 = load("cyp_knocks/cyp_33000.0_1.jld")
t = data1["t"]
u1 = data1["u"]
t1 = data1["t"]
u2 = data2["u"]
t2 = data2["t"]
u3 = data3["u"]
t3 = data3["t"]

####################### Generate Plots

## Boundary Sharpening Diagrams

gₕ = u1[2][:,31:61,7]
gₖ = u1[2][:,31:61,8]

PlotUtils.default_cgrad(default = :sequential, sequential = :redsblues)

p11 = heatmap(x,y,gₕ,cbar =  false,
              ylabel="Wildtype",clims=(0.0,2.0))

              p11 = heatmap(x,y,gₕ,cbar =  true,
                            ylabel="Wildtype",clims=(0.0,2.0))

gₕ = u1[19][:,31:61,7] # 1/2 hour
gₖ = u1[19][:,31:61,8]

p12 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u1[37][:,31:61,7] # 1 hour
gₖ = u1[37][:,31:61,8]

p13 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u1[2*37-1][:,31:61,7] # 2 hours
gₖ = u1[2*37-1][:,31:61,8]

p14 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

p1 = plot(p11,p12,p13,p14,layout=grid(1,4))

gₕ = u2[2][:,31:61,7]
gₖ = u2[2][:,31:61,8]

p21 = heatmap(x,y,gₕ,cbar =  false, ylabel="BP Knockdown",clims=(0.0,2.0))

gₕ = u2[19][:,31:61,7] # 1/2 hour
gₖ = u2[19][:,31:61,8]

p22 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u2[37][:,31:61,7] # 1 hour
gₖ = u2[37][:,31:61,8]

p23 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u2[2*37-1][:,31:61,7] # 2 hours
gₖ = u2[2*37-1][:,31:61,8]

p24 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u3[2][:,31:61,7]
gₖ = u3[2][:,31:61,8]

p31 = heatmap(x,y,gₕ,cbar =  false, ylabel="Cyp Knockdown",clims=(0.0,2.0))

gₕ = u3[19][:,31:61,7] # 1/2 hour
gₖ = u3[19][:,31:61,8]

p32 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u3[37][:,31:61,7] # 1 hour
gₖ = u3[37][:,31:61,8]

p33 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

gₕ = u3[2*37-1][:,31:61,7] # 2 hours
gₖ = u3[2*37-1][:,31:61,8]

p34 = heatmap(x,y,gₕ,cbar =  false,clims=(0.0,2.0))

p = plot(p11,p12,p13,p14,p21,p22,p23,p24,p31,p32,p33,p34,layout=grid(3,4))

####################### Calculate boundary variances

#### Frequency

boundary_stds_bp100 = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("wildtype/bp_100.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds4(u)
  boundary_stds_bp100[:,i] = boundary_stds
end
p2 = plot(t[1:73]/60/60 .+ 10,mean(boundary_stds_bp100,2)[1:73],lw=3,
        label="Wildtype",xlabel="Time (hpf)",
        ylabel="# Displaced\nCells",
        yaxis = [0,Inf])
boundary_stds_bp5 = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("bp_knocks/bp_5.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds4(u)
  boundary_stds_bp5[:,i] = boundary_stds
end
plot!(p2,t[1:73]/60/60 .+ 10,mean(boundary_stds_bp5,2)[1:73],
      lw=3,label="BP Knockdown")

boundary_stds_cyp = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("cyp_knocks/cyp_33000.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds4(u)
  boundary_stds_cyp[:,i] = boundary_stds
end
plot!(p2,t[1:73]/60/60 .+ 10,mean(boundary_stds_cyp,2)[1:73],
      lw=3,label="Cyp Knockdown")

### Maximum Length

boundary_stds_bp100 = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("wildtype/bp_100.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds3(u,maximum)
  boundary_stds_bp100[:,i] = boundary_stds
end

p3 = plot(t[1:73]/60/60 .+ 10,mean(boundary_stds_bp100,2)[1:73]./10,
          lw=3,legend=false,
          xlabel="Time (hpf)", ylabel="Max Boundary\nLength\n(cell lengths)",
          yaxis = [0,12])

boundary_stds_bp5 = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("bp_knocks/bp_5.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds3(u,maximum)
  boundary_stds_bp5[:,i] = boundary_stds
end
plot!(p3,t[1:73]/60/60 .+ 10,mean(boundary_stds_bp5,2)[1:73]./10,lw=3)

boundary_stds_cyp = Matrix{Float64}(length(t),10)
for i in 1:10
  data = load("cyp_knocks/cyp_33000.0_$(i).jld"); u = data["u"]
  t = data["t"]
  boundary_stds = calculate_boundary_stds3(u,maximum)
  boundary_stds_cyp[:,i] = boundary_stds
end
plot!(p3,t[1:73]/60/60 .+ 10,mean(boundary_stds_cyp,2)[1:73]./10,lw=3)

l = @layout [  a{0.75h}
             [b c]]

plot(p,p2,p3,layout=l,size=(750,500))

################################################################################

## Noise Success

Vbps = 5.0:10.0:105.0
ϵₕₖs = 0.1:0.025:0.3
num_sims = length(Vbps)*length(ϵₕₖs)
σRAs = Vector{Float64}(num_sims)
σGenes = Vector{Float64}(num_sims)
boundary_sharpness = Vector{Float64}(num_sims)
mean_boundary_sharpness = Vector{Float64}(num_sims)
number_boundary_over = Vector{Float64}(num_sims)
iter = 0
t[73]./60./60 == 2 #hpf
for Vbp in Vbps, ϵₕₖ in ϵₕₖs
  data = load("noise_noise_plot/bp_$(Int(Vbp))_epshk_$(ϵₕₖ).jld")
  u = data["u"]
  iter += 1
  σRAs[iter] = effective_σRA(u)
  boundary_sharpness[iter] = calculate_boundary_stds3(u,maximum)[73]
  mean_boundary_sharpness[iter] = calculate_boundary_stds3(u,mean)[73]
  number_boundary_over[iter] = calculate_boundary_stds4(u)[73]
  σGenes[iter] = ϵₕₖ
end

sharpened = boundary_sharpness.< 30

p1 = scatter(σGenes[sharpened],σRAs[sharpened],
        yscale=:log10,markersize=10,markercolor=:green,
        title="Maximum Displacement < 3 Cell Lengths",
        xlabel="Gene Regulatory Noise",
        ylabel="Effective RA Signal Noise",
        labels="Successful Shapening Event")
scatter!(p1,σGenes[!sharpened],σRAs[!sharpened],yscale=:log10,markersize=10,
        yscale=:log10,markersize=10,markercolor=:red,
        labels="Unsuccessful Shapening Event")

mean_sharpened = mean_boundary_sharpness .< 5

p2 = scatter(σGenes[mean_sharpened],σRAs[mean_sharpened],
        yscale=:log10,markersize=10,markercolor=:green,
        title="Mean Displacement < 1/2 Cell Lengths",
        xlabel="Gene Regulatory Noise",
        ylabel="Effective RA Signal Noise",
        labels="Successful Shapening Event")
scatter!(p2,σGenes[!mean_sharpened],σRAs[!mean_sharpened],
        yscale=:log10,markersize=10,markercolor=:red,
        labels="Unsuccessful Shapening Event")

number_sharpened = number_boundary_over .<= 3

plot(p1,p2)

p5 = scatter(σGenes[number_sharpened],σRAs[number_sharpened],
        yscale=:log10,markersize=10,markercolor=:green,
        xlabel="Gene Regulatory Noise",
        ylabel="Effective RA Signal Noise",
        labels="Successful Shapening Event")
scatter!(σGenes[!number_sharpened],σRAs[!number_sharpened],
         yscale=:log10,markersize=10,markercolor=:red,
         labels="Unsuccessful Shapening Event")

################################################################################

## Cyp Shift

RAR_wt = u1[end][:,31:61,4]
RAR_cyp = u3[end][:,31:61,4]
plot(X[1,31:61],mean(RAR_wt,1)',lw=3,labels="Wildtype",xlabel="Anterior-Posterior (um)",ylabel="Average [RAR] (uM)")
plot!(X[1,31:61],mean(RAR_cyp,1)',lw=3,labels="Cyp Knockdown")

################################################################################

## Supplement: Sharpening Index

boundary_stds_bp100 = Matrix{Float64}(length(t),10)
for i in 1:10
 data = load("wildtype/bp_100.0_$(i).jld"); u = data["u"]
 t = data["t"]
 boundary_stds = calculate_boundary_stds2(u)
 boundary_stds_bp100[:,i] = boundary_stds
end
plot(t,mean(boundary_stds_bp100,2),lw=3,label="Wildtype",xlabel="Time", ylabel="Sharpness Index")
boundary_stds_bp5 = Matrix{Float64}(length(t),10)
for i in 1:10
 data = load("bp_knocks/bp_5.0_$(i).jld"); u = data["u"]
 t = data["t"]
 boundary_stds = calculate_boundary_stds2(u)
 boundary_stds_bp5[:,i] = boundary_stds
end
plot!(t,mean(boundary_stds_bp5,2),lw=3,label="BP Knockdown")
boundary_stds_cyp = Matrix{Float64}(length(t),10)
for i in 1:10
 data = load("cyp_knocks/cyp_33000.0_$(i).jld"); u = data["u"]
 t = data["t"]
 boundary_stds = calculate_boundary_stds2(u)
 boundary_stds_cyp[:,i] = boundary_stds
end
plot!(t,mean(boundary_stds_cyp,2),lw=3,label="Cyp Knockdown")

####################### Gradient Plots

p1 = surface(X[:,31:61],Y[:,31:61],u1[end][:,31:61,2],title="[RAin]")
p2 = surface(X[:,31:61],Y[:,31:61],u1[end][:,31:61,4],zlims=[0.0,2.1],clims=(0.0,2.1),title="RAR")
p3 = surface(X[:,31:61],Y[:,31:61],u2[end][:,31:61,2],title="[RAin], BP Knockdown")
p4 = surface(X[:,31:61],Y[:,31:61],u2[end][:,31:61,4],zlims=[0.0,2.1],title="[RAR], BP Knockdown")

plot(p1,p2,p3,p4,cbar=false)

###################### Gradient Heatmaps

p1 = heatmap(x,y,u1[end][:,31:61,2],title="RAin",cbar=false)
p2 = heatmap(x,y,u1[end][:,31:61,4],zlims=[0.0,2.1],cbar=false,clims=(0.0,2.1),title="RAR")

p3 = heatmap(x,y,u2[end][:,31:61,2],cbar=false,title="RAin")
p4 = heatmap(x,y,u2[end][:,31:61,4],cbar=false,zlims=[0.0,2.1],title="RAR")

plot(p1,p2,p3,p4,cbar=false)
