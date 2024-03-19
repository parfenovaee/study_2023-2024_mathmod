using Plots
using DifferentialEquations

M0_1 = 2.5
M0_2 = 1.8
p_c = 20.0
N = 23.0
q = 1.0
tau_1 = 16.0
tau_2 = 19.0
p_1 = 13.0
p_2 = 11.0

a_1 = p_c/(tau_1*tau_1*p_1*p_1*N*q)
a_2 = p_c/(tau_2*tau_2*p_2*p_2*N*q)
b = p_c/(tau_1*tau_1*p_1*p_1*tau_2*tau_2*p_2*p_2*N*q)
c_1 = (p_c-p_1)/(tau_1*p_1)
c_2 = (p_c-p_2)/(tau_2*p_2)

start = [M0_1, M0_2]
timee = (0.0, 30.0)

function one_fun(du, u, p, t)
   du[1] = u[1] - b/c_1*u[1]*u[2]-a_1/c_1*u[1]*u[1]
   du[2] = c_2/c_1*u[2] - b/c_1*u[1]*u[2]-a_2/c_1*u[2]*u[2]
end

equat = ODEProblem(one_fun, start, timee)
solv = solve(equat, dtmax=0.01)

M_1 = [u[1] for u in solv.u]
M_2 = [u[2] for u in solv.u]

plot1 = plot(dpi = 600, legend =:bottomright, bg =:white, title="Изменение оборотных средств в 1 случае")
plot!(plot1, solv.t, M_1, label="Изменения объемов продаж 1 фирмы", color =:green)
plot!(plot1, solv.t, M_2, label="Изменения объемов продаж 2 фирмы", color =:blue)

savefig(plot1, "lab08_1.png")
