include(joinpath(@__DIR__, "config.jl"))

# Fix... x1 doesn't participate...

function solve_model_bi(opt_factory, run_name)
    println("Complete solving problem ($(run_name)) with arity = 2")
    m = Model(opt_factory)
    @variable(m,  0.1 <= x <= 1.0)
    @NLobjective(m, Min, (-(x - x^2))*(-(x^3 - exp(x))))
    optimize!(m)
    return m
end

function solve_model_tri_p(opt_factory, run_name)
    println("Complete solving problem ($(run_name)) with arity = 3")
    m = Model(opt_factory)
    @variable(m,  -0.5 <= x <= 1.0)
    @NLobjective(m, Min, 1.0*((x - x^2)*(x^3 - exp(x))))
    optimize!(m)
    return m
end

function solve_model_tri(opt_factory, run_name)
    println("Complete solving problem ($(run_name)) with arity = 3 routine")
    m = Model(opt_factory)
    @variable(m,  -0.5 <= x <= 1.0)
    @NLobjective(m, Min, 1.0*(x - x^2)*(x^3 - exp(x)))
    optimize!(m)
    return m
end

function f(xp, xpr, xL, xU)

    xr = MC{1,NS}(xpr,Interval(xL,xU),1)
    t1r = (xr - xr^2)
    t2r = (xr^3 - exp(xr))
    t3r = t1r*t2r

    x = MC{1,NS}(xp,Interval(xL,xU),1)
    t1 = (x - x^2)
    t2 = (x^3 - exp(x))
    t3 = t1*t2

    println("t1 = $t1, t1r = $(t1r)")
    println("t2 = $t2, t3r = $(t2r)")
    println("t3 = $t3, t3r = $(t3r)")
    return 
end

println("Running Benchmark - Test Function 2")
#m_eago        = solve_model_bi(eago, :eago)
#m_eago        = solve_model_tri(eago, :eago)
#m_eago        = solve_model_tri_p(eago, :eago)

#m_eago_grad   = solve_model_bi(eago_grad, :eago_grad)
#m_eago_grad   = solve_model_tri(eago_grad, :eago_grad)
#m_eago_grad   = solve_model_tri_p(eago_grad, :eago_grad)

#m_eago_enum   = solve_model_bi(eago_enum, :eago_enum)
#m_eago_enum   = solve_model_tri(eago_enum, :eago_enum)
#m_eago_enum   = solve_model_tri_p(eago_enum, :eago_enum)

m_eago_affine = solve_model_bi(eago_affine, :eago_affine)


#m_scip        = solve_model(scip, :scip)

#f(, , 0.8428052734375, 1.0)