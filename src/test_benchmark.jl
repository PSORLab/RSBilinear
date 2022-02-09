include(joinpath(@__DIR__, "config.jl"))

function generate_model(nvar, v)

    Q = zeros(3*nvar, 3*nvar)
    for i in 1:3*nvar, j in 1:3*nvar
        if i <= j
            if rand() <= v
                Q[i,j] = 2*rand() - 1
            end
        end
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
                     m = Model()\n
                     @variable(m, -1 <= x[i=1:$nvar] <= 1)\n
                     @variable(m, $(obj_bnd.lo) <= q <= $(obj_bnd.hi))\n
                     add_NL_constraint(m, :("*obj*" - \$q <= 0.0))\n
                     @objective(m, Min, q)\n
                     return m\n
                    "
              )
end

function create_lib()
end