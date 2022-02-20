include(joinpath(@__DIR__, "config.jl"))

# Fix... x1 doesn't participate...

function solve_model(opt_factory, run_name)
    m = Model(opt_factory)
    @variable(m,  0.1 <= x <= 0.9)
    @variable(m,  0.1 <= y <= 0.9)
    @NLobjective(m, Min, 1.0*(x^2 - x)*(y^2 - y))
    optimize!(m)
    println("Complete solving problem ($(run_name)) in get")
    return m
end

println("Running Benchmark - Test Function 1")
m_eago        = solve_model(eago, :eago)
m_eago_grad   = solve_model(eago_grad, :eago_grad)
#m_eago_enum   = solve_model(eago_enum, :eago_enum)
#m_eago_affine = solve_model(eago_affine, :eago_affine)
#m_scip        = solve_model(scip, :scip)