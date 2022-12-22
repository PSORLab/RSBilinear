include(joinpath(@__DIR__, "config.jl"))

function solve_model(opt_factory, run_name)
    m = Model(opt_factory)

    # i = 1 is placeholder to simplify indexing
    xL = [nothing; 7.0; 0.1; -1.078; -1.078]
    xU = [nothing; ; 0.3; 1.078; 1.078]
    @variable(m, xL[i-1] <= x[i=2:5] <= xU[i-1])
    @variable(m, z[i=1:2], Bin)

    @constraint(m, z[1] + z[2] <= 1)

    c = [13.1; 0.021]
    v = 5.0

    function objective_interval_bound(x, z)
        y1ex = 3.55 + 0.27*c[1] + 0.58*c[2] + 60.6*x[2] - 2.8*c[1]*x[2] - 2.3*c[2]*x[2]
        y2ex = 126585.5 - 21466*y1ex + 520.43*x[3] + 56.29*z[1] + 315.95*z[2] - 43.72*x[3]*y1ex + 3.74*x[3]^2 + 
               910.1*y1ex^2 - 30.6*y1ex*z[1] - 173.17*y1ex*z[2]
        y3ex = 9.16 + 0.092*x[3] + 0.73*x[4] + 0.64*x[3]*x[4] - 0.49*x[4]^2 - 0.13*x[4]*y2ex + 0.0019*y2ex^2 + 0.018*y2ex
    
        return (y3ex_1 - v)^2
    end

    xIntv = Interval.(xL, xU)
    zIntv = [Interval(0,1) for i=1:2]
    fIntv = objective_interval_bound(xIntv, zIntv)
    @variable(m, fIntv.lo <= η <= fIntv.hi)

    y1ex = :(3.55 + 0.27*$(c[1]) + 0.58*$(c[2]) + 60.6*$(x[2]) - 2.8*$(c[1])*$(x[2]) - 2.3*$(c[2])*$(x[2]))
    y2ex = :(126585.5 - 21466*$y1ex + 520.43*$(x[3]) + 56.29*$(z[1]) + 315.95*$(z[2]) - 43.72*$(x[3])*$(y1ex) + 3.74*$(x[3])^2 + 
           910.1*$y1ex^2 - 30.6*$y1ex*$(z[1]) - 173.17*$y1ex*$(z[2]))
    y3ex = :(9.16 + 0.092*$(x[3]) + 0.73*$(x[4]) + 0.64*$(x[3])*$(x[4]) - 0.49*$(x[4])^2 - 0.13*$(x[4])*$y2ex + 0.0019*$y2ex^2 + 0.018*$y2ex)
    add_NL_constraint(m, :(($y3ex_1 - $v)^2 - $η <= 0.0))
    @objective(m, Min, η)

    optimize!(m)
    println("Complete solving problem ($(run_name)) in get")
    println("Solve Time = $(solve_time(m))")
    return m
end

println("Running Benchmark - Response Surface Model")
#m_eago        = solve_model(eago, :eago)
m_eago_grad   = solve_model(eago_grad, :eago_grad)
#m_eago_enum   = solve_model(eago_enum, :eago_enum)
#m_eago_affine = solve_model(eago_affine, :eago_affine)
#m_scip        = solve_model(scip, :scip)
