u01  = GpuArray(zeros(M,N))
u02  = GPUArray(zeros(M,N))
u03  = GPUArray(0.2ones(M,N))
u04  = GPUArray(zeros(M,N))

u0 = ArrayPartition(u01,u02,u03,u04)
function hox_krox(du,u,p,t)
  RAout,RAin,gₕ,gₖ = u.x
  dRAout,dRAin,dgₕ,dgₖ = du.x
  Cyp .= !NCYP_kmax.*(kdeg.*γ.*RAin.^n)./(1.+γ.*RAin.^n.+f₀.*exp.(-λ.*(xf.-X))) .+
         NCYP_kmax.*kmax
  dRAout .= DRA.*(RAout*Ax + Ay*RAout) .+ VRA .- (1.+β).*kA.*RAout .+ kA.*RAin
  dRAin  .= kA.*RAout .- (kA .+ Cyp).*RAin
  if t > 500
    dgₕ .= (cₕ.*gₕ.^nₕ .+ (κₕ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₕ*RAin).^m) - dₕ.*gₕ
    dgₖ .= (cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) ./ (1 .+ cₕ.*gₕ.^nₕ .+ cₖ.*gₖ.^nₖ .+ (κₖ*RAin).^m) - dₖ.*gₖ
  end
  nothing
end
function hox_krox_noise(du,u,p,t)
  RAout,RAin,gₕ,gₖ = u.x
  dRAout,dRAin,dgₕ,dgₖ = du.x
  dRAout .= ϵout .* RAout
  dRAin  .= ϵin  .* RAin
  if t > 500
    dgₕ .= aₕ .* gₕ
    dgₖ .= aₖ .* gₖ
  end
  nothing
end

prob = SDEProblem(hox_krox,hox_krox_noise,u0,(0.0,1000.0))
sol = solve(prob,SRIW1())
