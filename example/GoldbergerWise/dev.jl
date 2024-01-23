using ModelingToolkit, DifferentialEquations
using IntervalArithmetic, IntervalRootFinding

@parameters a
@variables t, x(t)
D = Differential(t)
eqs = [D(x) ~ a * x]
@named sys = ODESystem(eqs)
sys = structural_simplify(sys)
tspan = (0..0,2..2)
function f(p0)
    u0 = [x => 1.0..1.0]
    p = [a => p0]
    prob = ODEProblem(sys, u0, tspan, p)
    return solve(prob, Tsit5(), dt=1e-2)(0)[1]
end
roots(f, 0..1)

x = complex(0,-pi..pi)
sinh(x)

@parameters σ ρ β
@variables t x(t) y(t) z(t)
D = Differential(t)

eqs = [D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z]

@named sys = ODESystem(eqs)
sys = structural_simplify(sys)
u0 = [D(x) => 2.0,
    x => 1.0..2.0,
    y => 0.0,
    z => 0.0]

p = [σ => 28.0,
    ρ => 10.0,
    β => 8 / 3]

tspan = (0.0..0.0, 10.0..10.0)
prob = ODEProblem(sys, u0, tspan, p, jac = true)
sol = solve(prob, ImplicitEuler())
using Plots;
plot(sol, idxs = (x, y))