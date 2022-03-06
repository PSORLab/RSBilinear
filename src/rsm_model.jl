include(joinpath(@__DIR__, "config.jl"))

# Fix... x1 doesn't participate...

function solve_model(opt_factory, run_name)
    m = Model(opt_factory)

    xL = [1200.0; 7.0; 0.1; -1.078; -1.078]
    xU = [1500.0; 10.0; 0.3; 1.078; 1.078]
    @variable(m, xL[i] <= x[i=2:5] <= xU[i])
    @variable(m, z[i=1:2], Bin)

    @constraint(m, z[1] + z[2] <= 1)

    c = [13.1; 0.021]

    ct1 = -0.237*c[1] + 0.189*c[2] 
    ct2 = -0.359*c[1] - 0.0004*c[2]
    ct3 = -0.047*c[1] + 0.122*c[2]

    coeff1 = 3.283 - 0.359*c[1] - 0.0004*c[2]
    coeff2 = -3.01*10^(-6)
    coeff3 = 0.0073 - 1.007*c[2]
    v = [12, 30]

    function objective_interval_bound(x, z)
        y1ex_1 = 14.84 + coeff1*x[3] + ct1
        #y1ex_2 = 0.576 + coeff2*x[1] + coeff3*x[3] + ct3
        y1ex_3 = -0.006 + 0.006*x[2] + 0.0038*x[3] - 0.0008*x[2]*x[3]
    
        y2ex_1 = 0.9015 + 0.4874*x[4] + 5.3102*y1ex_1 + 0.0067*z[1] - 0.03402*z[2] + 0.0438*x[4]*z[1] + 0.276*x[4]*z[2] - 0.67*x[4]^2 - 0.39*y1ex_1^2
        y2ex_2 = -0.2 - 0.159*x[4] - 1.183*y1ex_1 - 0.1317*x[4]^2 + 0.1057*y1ex_1^2
    
        y3ex_1 = 0.034 + 0.00137*y2ex_1 + 0.0617*x[5] - 0.0077*y2ex_1*x[5] + 0.0103*x[5]^2
        y3ex_2 = -12.77 + 2.457*y2ex_1 - 0.927*y2ex_2 + 2.876*x[5] + 2.853*x[5]^2
    
        return (y3ex_1 - v[1])^2 + (y1ex_3 + y2ex_2 + y3ex_2 - v[2])^2
    end

    xIntv = Interval.(xL, xU)
    zIntv = [Interval(0,1) for i=1:2]
    fIntv = objective_interval_bound(xIntv, zIntv)
    @variable(m, fIntv.lo <= η <= fIntv.hi)

    y1ex_1 = :(14.84 + $coeff1*($(x[3])) + $ct1)
    #y1ex_2 = :(0.576 + $coeff2*($(x[1])) + $coeff3*($(x[3])) + $ct3)
    y1ex_3 = :(-0.006 + 0.006*($(x[2])) + 0.0038*($(x[3])) - 0.0008*($(x[2]))*($(x[3])))
    
    y2ex_1 = :(0.9015 + 0.4874*($(x[4])) + 5.3102*($y1ex_1) + 0.0067*($(z[1])) - 0.03402*($(z[2])) + 0.0438*(($(x[4]))*($(z[1]))) + 0.276*(($(x[4]))*($(z[2]))) - 0.67*($(x[4]))^2 - 0.39*($y1ex_1)^2)
    y2ex_2 = :(-0.2 - 0.159*($(x[4])) - 1.183*($y1ex_1) - 0.1317*($(x[4]))^2 + 0.1057*($y1ex_1)^2)
    
    y3ex_1 = :(0.034 + 0.00137*($y2ex_1) + 0.0617*($(x[5])) - 0.0077*(($y2ex_1)*($(x[5]))) + 0.0103*($(x[5]))^2)
    y3ex_2 = :(-12.77 + 2.457*($y2ex_1) - 0.927*($y2ex_2) + 2.876*($(x[5])) + 2.853*($(x[5]))^2)
    
    add_NL_constraint(m, :(($y3ex_1 - $(v[1]))^2 + ($y1ex_3 + $y2ex_2 + $y3ex_2 - $(v[2]))^2 - $η <= 0.0))
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
