using Plots, EAGO

X = Interval(0.5, 0.8625)
Y = Interval(-1.0, 0.45)

#X = Interval(0.5, 0.9)
#Y = Interval(-1.0, 0.5)

#X = Interval(-1.0, 0.5)
#Y = Interval(-1.0, 1.0)

x0 = mid(X)
y0 = mid(Y)

xrng = X.lo:0.02:X.hi
yrng = Y.lo:0.02:Y.hi

m = 2
n = 3

# NEXT: ADD AFFINE CUTS TO INTERMEDIATE... (Not needed...)
# NEXT: CHECK SUBGRADIENTS AND REGENERATE CUTS...

#=
function estimator_extrema()
    x0mc = MC{2,NS}(x0, X, 1)
    y0mc = MC{2,NS}(y0, Y, 2)
    t1a = x0mc^m - x0mc
    t2a = y0mc^n - y0mc
    u1 = t1a.cv + t1a.cv_grad[1]*(X - x0) + t1a.cv_grad[2]*(Y - y0)
    u2 = t2a.cv + t2a.cv_grad[1]*(X - x0) + t2a.cv_grad[2]*(Y - y0)
    v1 = t1a.cc + t1a.cc_grad[1]*(X - x0) + t1a.cc_grad[2]*(Y - y0)
    v2 = t2a.cc + t2a.cc_grad[1]*(X - x0) + t2a.cc_grad[2]*(Y - y0)
    v1n = -v1
    v2n = -v2
    u1.hi, u2.hi, v1n.hi, v2n.hi
end

function estimator_under(t1, t2, x, y)
    x0mc = MC{2,NS}(x0, X, 1)
    y0mc = MC{2,NS}(y0, Y, 2)
    t1a = x0mc^m - x0mc
    t2a = y0mc^n - y0mc
    u1 = t1a.cv + t1a.cv_grad[1]*(x - x0) + t1a.cv_grad[2]*(y - y0)
    u2 = t2a.cv + t2a.cv_grad[1]*(x - x0) + t2a.cv_grad[2]*(y - y0)
    u1, u2, t1a.cv_grad, t2a.cv_grad
end

function estimator_over(t1, t2, x, y)
    x0mc = MC{2,NS}(x0, X, 1)
    y0mc = MC{2,NS}(y0, Y, 2)
    t1 = x0mc^m - x0mc
    t2 = y0mc^n - y0mc
    v1 = t1.cc + t1.cc_grad[1]*(x - x0) + t1.cc_grad[2]*(y - y0)
    v2 = t2.cc + t2.cc_grad[1]*(x - x0) + t2.cc_grad[2]*(y - y0)
    -v1, -v2, -t1.cc_grad, -t2.cc_grad
end
=#

#=
function estimator_extrema()
    x0l_mc = MC{2,NS}(X.lo, X, 1)
    x0u_mc = MC{2,NS}(X.hi, X, 1)
    y0l_mc = MC{2,NS}(Y.lo, Y, 2)
    y0u_mc = MC{2,NS}(Y.hi, Y, 2)
    t1_l = x0l_mc^m - x0l_mc
    t1_u = x0u_mc^m - x0u_mc
    t2_l = y0l_mc^n - y0l_mc
    t2_u = y0u_mc^n - y0u_mc
    u1max = max(t1_l.cv, t1_u.cv)
    u2max = max(t2_l.cv, t2_u.cv)
    v1nmax = min(t1_l.cc, t1_u.cc)
    v2nmax = min(t2_l.cc, t2_u.cc)
    u1max, u2max, -v1nmax, -v2nmax
end

function estimator_under(t1, t2, x, y)
    t1.cv, t2.cv, t1.cv_grad, t2.cv_grad
end

function estimator_over(t1, t2, x, y)
    -t1.cc, -t2.cc, -t1.cc_grad, -t2.cc_grad
end
=#

function estimator_extrema()
    x0aff = EAGO.AffineEAGO{2}(x0, X, 1)
    y0aff = EAGO.AffineEAGO{2}(y0, Y, 2)

    t1a = x0aff^m - x0aff
    t2a = y0aff^n - y0aff

    t1a_cv = t1a.c - t1a.Δ + 2.0*t1a.γ[1]*(X - x0)/(X.hi - X.lo) + 2.0*t1a.γ[2]*(Y - y0)/(Y.hi - Y.lo)
    t2a_cv = t2a.c - t2a.Δ + 2.0*t2a.γ[1]*(X - x0)/(X.hi - X.lo) + 2.0*t2a.γ[2]*(Y - y0)/(Y.hi - Y.lo)
    t1a_cc = -(t1a.c + t1a.Δ + 2.0*t1a.γ[1]*(X - x0)/(X.hi - X.lo) + 2.0*t1a.γ[2]*(Y - y0)/(Y.hi - Y.lo))
    t2a_cc = -(t2a.c + t2a.Δ + 2.0*t2a.γ[1]*(X - x0)/(X.hi - X.lo) + 2.0*t2a.γ[2]*(Y - y0)/(Y.hi - Y.lo))

    t1a_cv.hi, t2a_cv.hi, t1a_cc.hi, t2a_cc.hi
end

function estimator_under(t1, t2, x, y)
    x0aff = EAGO.AffineEAGO{2}(x0, X, 1)
    y0aff = EAGO.AffineEAGO{2}(y0, Y, 2)
    #@show x0aff
    #@show y0aff
    t1a = x0aff^m - x0aff
    t2a = y0aff^n - y0aff
 
    t1a_cv  = t1a.c - t1a.Δ + 2.0*t1a.γ[1]*(x - x0)/(X.hi - X.lo) + 2.0*t1a.γ[2]*(y - y0)/(Y.hi - Y.lo)
    t2a_cv  = t2a.c - t2a.Δ + 2.0*t2a.γ[1]*(x - x0)/(X.hi - X.lo) + 2.0*t2a.γ[2]*(y - y0)/(Y.hi - Y.lo)

    dZ = [X.hi - X.lo; Y.hi - Y.lo]

    t1a_cv_grad = SVector{2,Float64}(ntuple(i -> 2.0*t1a.γ[i]/dZ[i], Val(2)))
    t2a_cv_grad = SVector{2,Float64}(ntuple(i -> 2.0*t2a.γ[i]/dZ[i], Val(2)))
    t1a_cv, t2a_cv, t1a_cv_grad, t2a_cv_grad
end

function estimator_over(t1, t2, x, y)
    x0aff = EAGO.AffineEAGO{2}(x0, X, 1)
    y0aff = EAGO.AffineEAGO{2}(y0, Y, 2)
    t1a = x0aff^m - x0aff
    t2a = y0aff^n - y0aff
 
    t1a_cc  = t1a.c + t1a.Δ + 2.0*t1a.γ[1]*(x - x0)/(X.hi - X.lo) + 2.0*t1a.γ[2]*(y - y0)/(Y.hi - Y.lo)
    t2a_cc  = t2a.c + t2a.Δ + 2.0*t2a.γ[1]*(x - x0)/(X.hi - X.lo) + 2.0*t2a.γ[2]*(y - y0)/(Y.hi - Y.lo)

    dZ = [X.hi - X.lo; Y.hi - Y.lo]
    t1a_cc_grad = SVector{2,Float64}(ntuple(i -> -2.0*t1a.γ[i]/dZ[i], Val(2)))
    t2a_cc_grad = SVector{2,Float64}(ntuple(i -> -2.0*t2a.γ[i]/dZ[i], Val(2)))
    -t1a_cc, -t2a_cc, t1a_cc_grad, t2a_cc_grad
end

function f(x,y)
    xMC = MC{2,NS}(x, X, 1)
    yMC = MC{2,NS}(y, Y, 2)
    t1 = xMC^m - xMC
    t2 = yMC^n - yMC
    zIntv = t1.Intv*t2.Intv
    z = t1*t2
   # @show z
    u1max, u2max, v1nmax, v2nmax = estimator_extrema()
   @show u1max, u2max, v1nmax, v2nmax
    wIntv = z.Intv
    if (u1max < t1.Intv.hi) || (u2max < t2.Intv.hi)
        u1cv, u2cv, u1cvg, u2cvg = estimator_under(t1, t2, x, y)
        @show u1cv, u2cv, u1cvg, u2cvg 
        za_l = McCormick.mult_apriori_kernel(t1, t2, wIntv, u1cv, u2cv, u1max, u2max, u1cvg, u2cvg)
        #@show za_l
        z = z ∩ za_l
        #z = za_l
    end
    
    if (v1nmax > -t1.Intv.lo) || (v2nmax > -t2.Intv.lo)
        v1ccn, v2ccn, v1ccgn, v2ccgn = estimator_over(t1, t2, x, y)
        @show v1ccn, v2ccn, v1ccgn, v2ccgn 
        za_u = McCormick.mult_apriori_kernel(-t1, -t2, wIntv, v1ccn, v2ccn, v1nmax, v2nmax, v1ccgn, v2ccgn)
        @show za_u
        z = z ∩ za_u
    end
    @show z
    z
end

function fcv(x, y)
    zMC = f(x, y)
    zMC.cv
end

function fcc(x,y)
    zMC = f(x, y)
    zMC.cc
end

function fv(x,y)
    (x^m - x)*(y^n - y)
end

xr1 = 0.518125
yr1 = -0.9275
f0 = f(x0, y0)
@show f0
f1 = f(xr1, yr1)
function cv_g1(x, y)
    f0.cv + f0.cv_grad[1]*(x - x0) + f0.cv_grad[2]*(y - y0)
end
function cv_g2(x, y)
    f1.cv + f1.cv_grad[1]*(x - xr1) + f1.cv_grad[2]*(y - yr1)
end

function cc_g1(x, y)
    v = f0.cc + f0.cc_grad[1]*(x - x0) + f0.cc_grad[2]*(y - y0)
end
function cc_g2(x, y)
    v = f1.cc + f1.cc_grad[1]*(x - xr1) + f1.cc_grad[2]*(y - yr1)
end


surface(xrng, yrng, cv_g1)
surface!(xrng, yrng, cv_g2)
surface!(xrng, yrng, fcv)

#surface(xrng, yrng, fv)
#surface(xrng, yrng, fcc, camera = (20, 20))
#surface!(xrng, yrng, cc_g1, camera = (20, 20))
#surface!(xrng, yrng, cc_g2)

xlabel!("x")
ylabel!("y")
savefig("surface_stuff1")
println("ran surface gen")

println("ref 1")
f0 = f(x0, y0)
@show f0
println("ref 2")
f1 = f(xr1, yr1)

f0I = f0.cv + f0.cv_grad[1]*(xr1 - x0) + f0.cv_grad[2]*(yr1 - y0)
f1I = f0.cc + f0.cc_grad[1]*(xr1 - x0) + f0.cc_grad[2]*(yr1 - y0)
#f1I = f1.cv + f1.cv_grad[1]*(X - xr1) + f1.cv_grad[2]*(Y - yr1)
@show f0I
# MC{3, NS}(-0.5468342236328128, 0.5771286035156251, [-2.66832, 2.88813], [0.3955328125000001, -0.06036679687500004, -1.0], [0.3955328125000001, 0.15159375000000005, -1.0], false)
@show f1I
# MC{3, NS}(0.3092874172265627, 1.2738127189843755, [-2.66832, 2.88813], [-0.5256250000000001, 0.5194765625000001, -1.0], [0.3955328125000001, 0.15159375000000005, -1.0], false)
#ex1u, ex2u, ex1vn, ex2vn = estimator_extrema()

# IS AFFINE X GOOD?
# IS AFFINE Y GOOD?
# 

f1_test(x) = x^m - x
f2_test(x) = x^n - x

function estimator_under_x(x)
    x0aff = EAGO.AffineEAGO{1}(x0, X, 1)
    t1a = f1_test(x0aff)
    t1a_cv = t1a.c - t1a.Δ + 2*t1a.γ[1]*(x - mid(X))/(X.hi - X.lo)
    t1a_cv
end

function estimator_over_x(x)
    x0aff = EAGO.AffineEAGO{1}(x0, X, 1)
    t1a = f1_test(x0aff)
    t1a_cc  = t1a.c + t1a.Δ + 2*t1a.γ[1]*(x - mid(X))/(X.hi - X.lo)
    t1a_cc
end

function estimator_under_y(y)
    y0mc = EAGO.AffineEAGO{1}(y0, Y, 1)
    t2 = f2_test(y0mc)
    u2 = t2.c - t2.Δ + 2*t2.γ[1]*(y - mid(Y))/(Y.hi - Y.lo)
    u2
end

function estimator_over_y(y)
    y0mc = EAGO.AffineEAGO{1}(y0, Y, 1)
    t2 = f2_test(y0mc)
    v2 = t2.c + t2.Δ + 2*t2.γ[1]*(y - mid(Y))/(Y.hi - Y.lo)
    v2
end


function estimator_under_x_check(x, y)
    estimator_under(0.0, 0.0, x, y)
end


plt_x = plot(estimator_under_x, xrng)
plot!(plt_x, f1_test, xrng)
plot!(plt_x, estimator_over_x, xrng)

plt_y = plot(estimator_under_y, yrng)
plot!(plt_y, f2_test, yrng)
plot!(plt_y, estimator_over_y, yrng)

l = @layout [a ; c]

plot(plt_x, plt_y, layout = l)

savefig("subplot_fig")
