include(joinpath(@__DIR__, "config.jl"))

function populate_Q!(Q, nvar, v)
    Q = zeros(3*nvar, 3*nvar)
    for i in 1:3*nvar, j in 1:3*nvar
        if i <= j
            if rand() <= v
                Q[i,j] = 2*rand() - 1
            end
        end
    end
    return Q
end

function generate_model!(model_file_name, nvar, v)

    Q = zeros(3*nvar, 3*nvar)
    while iszero(sum(Q))
        populate_Q!(Q, nvar, v)
    end

    c = -rand(nvar).*510 .- 2

    @variables(m, -1 <= x[i=1:nvar] <= 1)
    y = Any[nothing i=1:3*nvar]
    ybnd = zeros(Interval{Float64}, 3*nvar)
    for i=1:nvar
        y[3*i-2] = "x[$i]^2"
        y[3*i-1] = "x[$i]^3"
        y[3*i] = "x[$i]^4"
        ybnd[3*i-2] = Interval(0,1)
        ybnd[3*i-2] = Interval(-1,1)
        ybnd[3*i-2] = Interval(0,1)
    end
    
    obj = ""
    obj_bnd = zero(Interval{Float64})
    for i in 1:3*nvar, j in 1:3*nvar
        if !iszero(Q[i,j])
            if first_nonzero
                obj = "$(Q[i,j])*("*y[i]*")*("*y[j]*")"
            else
                obj = obj*" + $(Q[i,j])*("*y[i]*")*("*y[j]*")"
            end
            obj_bnd += Q[i,j]*ybnd[i]*ybnd[j]
        end
    end
    for i=1:nvar
        obj = obj*" + $(c[i])*x[$i]"
    end
    obj_bnd += sum(c)

    open(model_file_name, "w") do file
        write(file, "using JuMP, EAGO\n
                     n = $nvar\n
                     v = $v\n
                     m = Model()\n
                     @variable(m, -1 <= x[i=1:$nvar] <= 1)\n
                     @variable(m, $(obj_bnd.lo) <= q <= $(obj_bnd.hi))\n
                     add_NL_constraint(m, :("*obj*" - \$q <= 0.0))\n
                     @objective(m, Min, q)\n
                     return m\n
                    "
              )
    println("wrote instance to $model_file_name")
end

function create_lib(file_path)
    nvec = [5, 7, 9]
    vvec = [0.3, 0.5, 0.7]
    for i = 1:200
        nvar = 
        v = rand(vvec)
        instance_file_path = file_path*"_$i_$nvar_$v"
        generate_model!(nvar, v)
    end
    return 
end

create_lib_path = joinpath(@__DIR__, "src", "MINLPLib.jl", "instances")

# set create_lib_files to true to generate new composite QPs in the problem libary
create_lib_files = false
create_lib_files && create_lib()

println("Begin running benchmark set.")
solvers = Dict{String,Any}()
solvers["SCIP--REG"] = scip
solvers["BARON-REG"] = baron

solvers["EAGO--REG"] = eago
solvers["EAGO--SUB"] = eago_grad
solvers["EAGO--AFF"] = eago_affine
solvers["EAGO--ENU"] = eago_enum

params = SolverBenchmarking.BenchmarkParams(time = 100, rerun = false, has_obj_bnd = false)
result_path = joinpath(@__DIR__, "benchmark_set_results")

SolverBenchmarking.run_solver_benchmark(result_path, solvers, "composite_qp", "composite_qp"; params = params)
SolverBenchmarking.summarize_results("composite_qp", result_path)

println("Finish running benchmark set.")