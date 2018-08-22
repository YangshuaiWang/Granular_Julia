at = gran.setup_cell(1000, 0.835, 2.0)
ϵ = [1.0 1.5; 1.5 0.5]*0.03
σ = [1.0 0.8; 0.8 0.88]*0.03
#ϵ = [1.0 1.0; 1.0 1.0]*0.05
#σ = [1.0 1.0; 1.0 1.0]*0.05

k = [1.0 1.0; 1.0 1.0]
α = [2.0 2.0; 2.0 2.0]

r1 = minimum(at.M)
r2 = maximum(at.M)

r0 = [2*r1 r1+r2; r1+r2 2*r2]

z = [1, 2]
#V = MultiLJ.MLJ(z, ϵ, σ)
V = MultiHZ.Multi_HZ(z, k, α, r0)
@show at.cell
@show energy(V, at)

set_calculator!(at, V)
myresult = minimise!(at, precond = :id, method = :lbfgs, gtol=1e-14, verbose=2)
