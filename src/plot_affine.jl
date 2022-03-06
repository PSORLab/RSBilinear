using EAGO, Plots

println(" ")
run_plots = true

l = -1.0
u = 1.0
x0 = l + 0.5*(u - l)
f(x) = 0.11664152193251676*((x^2 - x)*(x^3 - x)) #(x^2 - x)*(x^3 - x)


aff = f(EAGO.AffineEAGO{1}(x0, Interval(l,u), 1))
@show aff
@show f(Interval(l,u))

function aff_cv(x, x0, l, u)
    aff = f(EAGO.AffineEAGO{1}(x0, Interval(l,u), 1))
    aff.c - aff.Δ + aff.γ[1]*(x - x0)
end
aff_cv(x) = aff_cv(x, x0, l, u)

function aff_cc(x, x0, l, u)
    aff = f(EAGO.AffineEAGO{1}(x0, Interval(l,u), 1))
    aff.c + aff.Δ + aff.γ[1]*(x - x0)
end
aff_cc(x) = aff_cc(x, x0, l, u)

function mc_cv(x, x0, l, u)
    mc = f(MC{1,NS}(x0, Interval(l,u), 1))
    mc.cv + mc.cv_grad[1]*(x - x0)
end
mc_cv(x) = mc_cv(x, x0, l, u)

function mc_cc(x, x0, l, u)
    mc = f(MC{1,NS}(x0, Interval(l,u), 1))
    mc.cc + mc.cc_grad[1]*(x - x0)
end
mc_cc(x) = mc_cc(x, x0, l, u)

if run_plots
    xrng = l:0.01:u
    plot(f, xrng, label="f(x)",legend=:bottomright)
    plot!(aff_cv, xrng, label="affine(x) cv")
    plot!(aff_cc, xrng, label="affine(x) cc")
    #plot!(mc_cv, xrng, label="relax(x) cv")
    #plot!(mc_cc, xrng, label="relax(x) cc")
    savefig("affine_compare")
end
println("finished test file")