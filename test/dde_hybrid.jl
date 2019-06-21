using DiffEqJump, DelayDiffEq
using Test

doprint = true
# using Plots; pyplot()
doplot = false

# This example is due to @korsbo (#71) and has been modified as follows:
#   A VariableRateJump acts on u_1 and increases its value by 1
#   A ConstantRateJump acts on u_2 and decreases its value by 1

function dde(du, u, h, p, t)
  du[1] = h(p, t-p[1], idxs=1)
end

h(p, t; kwargs...) = 0
p = [ 5.]
tspan = (0., 10.)
u0 = [1., 100.]

prob = DDEProblem(dde, u0, h, tspan, p)

# add the jumps
prod_rate(u,p,t) = (1 + 9 * u[1] / (5 + u[1]))/10
prod_affect!(i) = (i.u[1] += 1)
prod_jump = VariableRateJump(prod_rate, prod_affect!)

deg_rate(u,p,t) = u[2]/10
deg_affect!(i) = (i.u[2] -= 1)
deg_jump = ConstantRateJump(deg_rate, deg_affect!)

### declare methods to check
# methods = (Direct(), FRM(), SortingDirect(), NRM(), RSSA(), DirectCR())
methods = (Direct(),)

if doplot
  plothand = plot(reuse=false, layout = grid(1,2), title = ["u1(t)" "u2(t)"], legend=nothing)
end

for method in methods
  jump_prob = JumpProblem(prob, Direct(), prod_jump, deg_jump)
  alg = MethodOfSteps(Tsit5())
  sol = solve(jump_prob, alg)

  if doplot
    plot!(plothand, sol)
  end

  if doprint
    println("   MethodOfSteps(Tsit5()) + ", typeof(method), ", sol[end] = ", sol[end])
  end

  @test sol[1,end] > u0[1]
  @test sol[2,end] < u0[2]
end
