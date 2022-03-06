
include(joinpath(@__DIR__, "config.jl"))

#=
Comments on formulation. Using x/y = x*y^(-1) reformulation as apriori relaxation
is implemented on * only currently and SCIP/baron don't support inv operator.
=#

n = 5
Δx = 1/(n-1)

# generate data for fitting the temperature data to (T_ref)
function compute_T(p::T) where T
    A  = zeros(T, n, n)
    b  = zeros(T, n)
    in_range_ref(i) = 31 <= i <= 65

    b[1] = 500
    for i = 2:(n-1)
        q0 = in_range_ref(i) ? 15000 : -3000
        b[i] = -(q0/p)*Δx^2
    end
    b[n] = 600

    q1 = -1.0
    for i = 2:(n - 1)
        A[i,i-1] = 1.0
        A[i,i+1] = 1.0
        A[i,i] = (-2.0 - (q1/p)*Δx^2)
    end
    A[1,1] = 1.0
    A[n,n] = 1.0
    A\b
end
p_ref = 0.5
T_ref = compute_T(p_ref)

@show T_ref[1]
@show T_ref[n]

# Maximum temperature specified in Mitsos 2009 to be 2000
obj = sum((T_ref[i] - Interval(0,2000))^2 for i=2:(n-1))

@show obj

# Min objective for SSE is always 0
obj_bnd = Interval(obj.lo, obj.hi)

#=
Adds constraints, objective, etc to model and solves
=#
function solve_model(opt_factory, run_name)
    m = Model(opt_factory)
    T1 = 500
    T101 = 600

    Δx = 1/(n-1)
    @variable(m, 0.01 <= p <= 10.0)
    
    in_range(i) = 51 <= i <= 61
    d_coeff = zeros(n)
    for i = 1:n
        q0 = in_range(i) ? 35000 : -500
        d_coeff[i] = -(isone(i) ? 500 : ((i == 101) ? 600 : q0))*Δx^2
    end
    #@NLexpression(m, d[i=1:n], d_coeff[i]*inv(p))
    d = Any[nothing for i=1:n]
    for i = 1:n
        d[i] =  add_NL_expression(m, :($(d_coeff[i])*$(p)^(-1)))
    end
    b = Any[nothing for i=1:n]
    b[1] = 1.0
    b[n] = 1.0
    for i = 2:(n-1)
        #b[i] = add_NL_expression(m, :((-2 + inv($(p))*$(Δx^2))))
        b[i] = add_NL_expression(m, :((-2 + $(Δx^2)*($(p))^(-1))))
    end

    cprime = Any[nothing for i=1:n]
    #cprime[1] = add_NL_expression(m, :(inv($(b[1]))))
    cprime[1] = add_NL_expression(m, :(($(b[1]))^(-1)))
    for i = 2:(n-1)
        #cprime[i] = add_NL_expression(m, :(inv($(b[i])) - $(cprime[i-1])))
        cprime[i] = add_NL_expression(m, :(($(b[i]) - $(cprime[i-1]))^(-1)))  
    end

    dprime = Any[nothing for i=1:n]
    #dprime[1] = add_NL_expression(m, :($(d[1])*inv($(b[1]))))
    dprime[1] = add_NL_expression(m, :($(d[1])*($(b[1]))^(-1)))
    for i = 2:n
        #dprime[i] = add_NL_expression(m, :(($(d[1]) - $(dprime[i-1]))*inv($(b[i]) - $(cprime[i-1]))))
        dprime[i] = add_NL_expression(m, :(($(d[i]) - $(dprime[i-1]))*($(b[i]) - $(cprime[i-1])^(-1))))
    end

    T = Any[nothing for i=1:n]
    T[n] = add_NL_expression(m, :($(dprime[n])))
    for i = (n - 1):-1:1
        T[i] = add_NL_expression(m, :($(dprime[i]) - $(cprime[i])*$(T[i + 1])))
    end
    @variable(m, obj_bnd.lo <= η <= obj_bnd.hi)
    @NLconstraint(m, sum((T[i] - T_ref[i])^2 for i=1:n) - η <= 0)
    @objective(m, Min, η)
    optimize!(m)
    println("Complete solving problem ($(run_name)) in get")

    return m
end

println("Running Benchmark - Heat Equation")
#m_eago        = solve_model(eago, :eago)
#m_eago_grad   = solve_model(eago_grad, :eago_grad)
#m_eago_enum   = solve_model(eago_enum, :eago_enum)
#m_eago_affine = solve_model(eago_affine, :eago_affine)
m_scip        = solve_model(scip, :scip)
#m_baron        = solve_model(baron, :baron)

function solve_model_nlexpr(opt_factory, run_name)
end
