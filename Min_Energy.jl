module Minimise_Energy

using JuLIP
using JuLIP.Potentials
using MultiHZ
using Pressure
import JuLIP: energy, forces, hessian_pos

export min_Energy

function min_Energy(at::AbstractAtoms)
    max_iter_num = 1000

    α = get_data(at, "Alpha")
    α = [α α; α α]
    k = 1.0
    k = [k k; k k]

    r = at.M
    r1 = minimum(r)
    r2 = maximum(r)
    r0 = [2*r1 r1+r2; r1+r2 2*r2]

    z = [1, 2]
    V = MultiHZ.Multi_HZ(z, k, α, r0)

    set_calculator!(at, V)

    myresult = minimise!(at, precond = :id, method = :lbfgs, gtol=1e-14, verbose=2)
    is_converged = myresult.f_converged || myresult.g_converged || myresult.x_converged
    optmsg = myresult

    iter_num = myresult.iterations
    f_calls = myresult.f_calls
    g_calls = myresult.g_calls 

    while  (myresult.iterations == max_iter_num) & !(is_converged)
        myresult = minimise!(at, precond = :id, method = :lbfgs, gtol=1e-14, verbose=2)
        iter_num += myresult.iterations
        f_calls += myresult.f_calls
        g_calls += myresult.g_calls 
    end
    optmsg = myresult
    optmsg.iterations = iter_num
    optmsg.f_calls = f_calls
    optmsg.g_calls = g_calls
    
    set_data!(at, "Optmsg", optmsg)
    
    Energy = energy(V,at)
    set_data!(at, "Energy", Energy)

    p = Pressure.cal_Pressure(at)
    set_data!(at, "Pressure", p)


    return at

end

end
