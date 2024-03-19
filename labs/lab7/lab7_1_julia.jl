using Plots
using DifferentialEquations

N = 810.0
n_0 = 11.0

function one(du, u, p, t)
   du[1] = (0.64 + 0.00014*u[1])*(N - u[1])
end

function two(du, u, p, t)
   du[1] = (0.000014 + 0.63*u[1])*(N - u[1])
end

function three(du, u, p, t)
   du[1] = (0.7*t + 0.4*cos(t)*u[1])*(N - u[1])
end

time = (0.0, 30.0)
timee = (0.0, 0.05)
start = [n_0]

equat1 = ODEProblem(one, start, time)
solv1 = solve(equat1, dtmax=0.01)  

equat2 = ODEProblem(two, start, timee)
solv2 = solve(equat2, dtmax=0.01)  

equat3 = ODEProblem(three, start, timee)
solv3 = solve(equat3, dtmax=0.01)  

n_1 = [u[1] for u in solv1.u]
n_2 = [u[1] for u in solv2.u]
n_3 = [u[1] for u in solv3.u]

max = 0;
max_t = 0;
max_n = 0;
for (i, t) in enumerate(solv2.t)
    if solv2(t, Val{1})[1] > max
        global max = solv2(t, Val{1})[1]
        global max_t = t
        global max_n = n_2[i]
    end
end

plot1 = plot(dpi = 300, legend= false, bg =:white, title="График распространения рекламы в 1 случае")
plot!(plot1, solv1.t, n_1, color =:green)

savefig(plot1, "lab07_1.png")

plot2 = plot(dpi = 300, legend= false, bg =:white, title="График распространения рекламы во 2 случае")
plot!(plot2, solv2.t, n_2, color =:green)
plot!(plot2, [max_t], [max_n], seriestype = :scatter, color = :green)

savefig(plot2, "lab07_2.png")

plot3 = plot(dpi = 300, legend= false, bg =:white, title="График распространения рекламы в 3 случае")
plot!(plot3, solv3.t, n_3, color =:green)

savefig(plot3, "lab07_3.png")