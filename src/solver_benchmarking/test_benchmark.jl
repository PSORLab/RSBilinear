# Adds packages that don't require special setup
using Pkg
#Pkg.develop(path = joinpath(@__DIR__, "McCormick.jl-master"))
Pkg.develop(path = joinpath(@__DIR__, "MINLPLib.jl"))
Pkg.develop(path = joinpath(@__DIR__, "EAGO.jl-SIPextension"))
Pkg.develop(path = joinpath(@__DIR__, "McCormick.jl-master"))
Pkg.add("JSON"); Pkg.add("DataFrames"); Pkg.add("CSV")
Pkg.add("JuMP"); Pkg.add("IntervalArithmetic"); Pkg.add("LaTeXStrings")
Pkg.add("Plots"); Pkg.add("StatsBase"); Pkg.add("BenchmarkProfiles")

# Loads solver_benchmarking module
include(joinpath(@__DIR__, "solver_benchmarking.jl"))
include(joinpath(@__DIR__, "config.jl"))

# Loads relevent modules
using JuMP, MINLPLib, SCIP, EAGO, IntervalArithmetic, GAMS
using CSV, DataFrames, LaTeXStrings, Plots, StatsBase, BenchmarkProfiles

function generate_model!(model_file_name, nvar, v)

    Q = zeros(3*nvar, 3*nvar)
    count = 0
    for i=1:1000
        for i = 1:(3*nvar), j = 1:(3*nvar)
            if i <= j
                t = rand()
                if t <= v
                    Q[i,j] = 2*rand() - 1
                end
            end
        end
        if !all(iszero, Q)
            break
        end
    end

    c = -rand(nvar).*510 .- 2

    y = Any[nothing for i=1:3*nvar]
    ybnd = zeros(Interval{Float64}, 3*nvar)
    for i=1:nvar
        y[3*i-2] = "\$(x[$i])^2"
        y[3*i-1] = "\$(x[$i])^3"
        y[3*i] = "\$(x[$i])^4"
        ybnd[3*i-2] = Interval(0,1)
        ybnd[3*i-1] = Interval(-1,1)
        ybnd[3*i] = Interval(0,1)
    end
    
    obj = ""
    obj_bnd = zero(Interval{Float64})
    first_nonzero = true
    for i in 1:3*nvar, j in 1:3*nvar
        if !iszero(Q[i,j])
            if first_nonzero
                obj = "$(Q[i,j])*(("*y[i]*")*("*y[j]*"))"
                first_nonzero = false
            else
                obj = obj*" + $(Q[i,j])*(("*y[i]*")*("*y[j]*"))"
            end
            obj_bnd += Q[i,j]*ybnd[i]*ybnd[j]
        end
    end
    for i=1:nvar
        obj = obj*" + $(c[i])*\$(x[$i])"
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
    end
    println("wrote instance to $model_file_name")
end

function create_lib(file_path)
    nvec = [10, 20, 30]
    vvec = [0.3, 0.5, 0.7]
    for i = 1:200
        nvar = rand(nvec)
        v = rand(vvec)
        instance_file_path = joinpath(file_path, "rs_quadratic_comp", "$(i)_$(nvar)_$(string(v)[1:2:3]).jl")
        println("generating model $instance_file_path")
        generate_model!(instance_file_path, nvar, v)
    end
    return 
end

# set create_lib_files to true to generate new ANNs in the problem libary
create_lib_files = false
create_lib_path = joinpath(@__DIR__, "MINLPLib.jl", "instances")
create_lib_files && create_lib(create_lib_path)


solvers = Dict{String,Any}()
solvers["SCIP--REG"] = scip
solvers["BARON-REG"] = baron
solvers["EAGO--REG"] = eago
#solvers["EAGO--SUB"] = eago_grad
#solvers["EAGO--AFF"] = eago_affine
#solvers["EAGO--ENU"] = eago_enum

function drop_num_underscore(x)
    y = replace(x, "_" => "")
    for i=0:9
        y = replace(y, "$i" => "")
    end
    return y
end

function print_summary_tables(df_comb, fig_name, env_lib, upper_limit)
    
    df_comb.SolveTime = map((x,y) -> ifelse(occursin("INFEASIBLE",x), upper_limit, y), df_comb.TerminationStatus, df_comb.SolveTime)
    df_comb.ShiftedSolveTime = df_comb.SolveTime .+ 1
    df_comb.CorrectlySolved = map(x -> (occursin("OPTIMAL",x) || occursin("LOCALLY_SOLVED",x)), df_comb.TerminationStatus)

    df_comb.IncorrectFeasibility = map(x -> occursin("INFEASIBLE",x), df_comb.TerminationStatus)
    # Presolve infeasibility in SCIP is labelled as OPTIMIZE_NOT_CALLED...
    df_comb.IncorrectFeasibility = map((x,y,z) -> ifelse(occursin("SCIP",x) & occursin("NOT_CALLED",y), true, z), df_comb.SolverName, df_comb.TerminationStatus, df_comb.IncorrectFeasibility)

    df_comb.ActFunc = map(drop_num_underscore, df_comb.InstanceName)

    gdf_comb_cs = groupby(df_comb, Symbol[:SolverName, :ActFunc])
    @show combine(gdf_comb_cs, :ShiftedSolveTime => x -> StatsBase.geomean(x) - 1)

    gdf_combt_correct = groupby(df_comb , Symbol[:SolverName, :CorrectlySolved, :ActFunc])
    @show combine(gdf_combt_correct, :ShiftedSolveTime => x -> count(fill(true,length(x))))

    df_comb_infeas = df_comb[df_comb.IncorrectFeasibility,:]
    gdf_check_infeasible  = groupby(df_comb_infeas, Symbol[:SolverName, :ActFunc])
    @show combine(gdf_check_infeasible, :SolveTime => x -> count(fill(true,length(x))))

    #gdf_status = groupby(df_comb, Symbol[:SolverName, :SolvedInTime])
    env_folder = joinpath(@__DIR__, "MINLPLib.jl", "instances", env_lib)
    name_anns = filter(x -> !occursin("gelu",x), readdir(joinpath(env_folder)))

    df_comb.UnsolvedRelativeGap = map((x,y,z) -> (~x & ~y & isfinite(z)), df_comb.CorrectlySolved, df_comb.IncorrectFeasibility, df_comb.RelativeGap)
    df_comb_gap = df_comb[df_comb.UnsolvedRelativeGap,:]
    gdf_comb_gap  = groupby(df_comb_gap, Symbol[:SolverName, :ActFunc])
    @show combine(gdf_comb_gap, :RelativeGap => x -> mean(x))

    trunc_solved_time = rand(200, 6)
    for (i,n) in enumerate(name_anns)
        plt_sdf = df_comb[occursin.(n[1:end-3], string.(df_comb.InstanceName)), :]
        trunc_solved_time[i,1] = plt_sdf[plt_sdf.SolveNamer .== "SCIP",:].SolveTime[1]  # SCIP
        trunc_solved_time[i,2] = plt_sdf[plt_sdf.SolveNamer .== "BARON",:].SolveTime[1] # BARON
        trunc_solved_time[i,3] = plt_sdf[plt_sdf.SolveNamer .== "EAGO",:].SolveTime[1]  # EAGO Env Entry
        trunc_solved_time[i,4] = plt_sdf[plt_sdf.SolveNamer .== "EAGO Sub",:].SolveTime[1]  # EAGO Exp Entry
        trunc_solved_time[i,5] = plt_sdf[plt_sdf.SolveNamer .== "EAGO Aff",:].SolveTime[1]  # EAGO Env Entry
        trunc_solved_time[i,6] = plt_sdf[plt_sdf.SolveNamer .== "EAGO Enum",:].SolveTime[1]  # EAGO Exp Entry
    end
    plt = performance_profile(PlotsBackend(), trunc_solved_time, ["SCIP", "BARON", "EAGO", "EAGO Sub", "EAGO Aff", "EAGO Enum"], linewidth = 1.5, linestyles=[:solid, :dash, :dashdot, :dot, :dashdot, :dot], legend=:bottomright)
    xlabel!("\$\\tau\$")
    xlims!(0.0,6.5)
    ylabel!("\$P(r_{p,s} \\leq \\tau : 1 \\leq s \\leq n_s)\$")
    ylims!(0.2,1.0)
    savefig(plt, joinpath(result_path, fig_name*".pdf"))
    show(plt)
end

result_path = joinpath(@__DIR__, "solver_benchmark_result")

params= SolverBenchmarking.BenchmarkParams(time = 300.0, rerun = false, has_obj_bnd = true)
SolverBenchmarking.run_solver_benchmark(result_path, solvers, "rs_quadratic_comp", "rs_quadratic_comp"; params = params)
SolverBenchmarking.summarize_results("rs_quadratic_comp", result_path)

df = DataFrame(CSV.File(joinpath(result_path, "rs_quadratic_comp", "result_summary.csv")))
print_summary_tables(df, "performance_profile", "rs_quadratic_comp", 300.0)
