using JuMP, EAGO

                     n = 10

                     v = 0.3

                     m = Model(EAGO.Optimizer)
                     set_optimizer_attribute(m, "output_iterations", 1)
                     set_optimizer_attribute(m, "verbosity", 4)
                     set_optimizer_attribute(m, "mul_relax_style", 1)
                     # 2 is affine arithmetic
                     # 3 is enumerated
                     set_optimizer_attribute(m, "subgrad_tighten", false)
                     set_optimizer_attribute(m, "iteration_limit", 200)

                    # FAILING APRIORI SG BOX... 0.68125, -0.275, 0.0 IN [0.5, -1.0, -2.0] TO [0.8625, 0.45, 2.0]

                     #@variable(m, 0.5 <= x <= 1)
                     #@variable(m, -1 <= y <= 1)
                     #@variable(m, -2.0 <= q <= 2.0)

                     @variable(m, 0.5 <= x <= 0.8625)
                     @variable(m, -1.0 <= y <= 0.45)
                     @variable(m, -2.0 <= q <= 2.0)

                     add_NL_constraint(m, :((($x^2 - $x)*($y^3 - $y)) - $q <= 0.0)) # FAILS ON EITHER... 2 or 3
                     #add_NL_constraint(m, :((($y^3 - $y) - $q <= 0.0))) # PASS
                     #add_NL_constraint(m, :((($x^2 - $x) - $q <= 0.0))) # PASS
                     #add_NL_constraint(m, :($x^3  - $q <= 0.0)) # PASS
                     #add_NL_constraint(m, :($x^2  - $q <= 0.0)) # PASSES?
                     #add_NL_constraint(m, :($x  - $q <= 0.0)) # PASSES?
                     @objective(m, Min, q)
                     optimize!(m)
                     @show objective_value(m)

function reference_relax(x0, l, u)
    x = [MC{3,NS}(x0[i], Interval(l[i], u[i]), i) for i=1:3]
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    x1_sqr = x[1]^2
    x2_cube = x[2]^3
    x1_sqr_m_x = x1_sqr - x1
    x2_cube_m_x = x2_cube - x2
    x1_times_x2 = x1_sqr_m_x*x2_cube_m_x

    println("x1 = $x1")
    println("x2 = $x2")
    println("x3 = $x3")

    println("x1_sqr = $x1_sqr in [$(l[1]^2), $(u[1]^2)]")
    println("x2_cube = $x2_cube in [$(min(l[2]^3,u[2]^3)), $(max(l[2]^3,u[2]^3))]")

    t3_min = minimum([-0.25, (l[1]^2-l[1]), (u[1]^2-u[1])])
    t3_max = maximum([(l[1]^2-l[1]), (u[1]^2-u[1])])
    println("x1_sqr_m_x = $x1_sqr_m_x in [$t3_min, $t3_max]")

    t4_min = minimum([-0.385, (l[1]^3-l[1]), (u[1]^3-u[1])])
    t4_max = maximum([0.385, (l[1]^3-l[1]), (u[1]^3-u[1])])
    println("x2_cube_m_x = $x2_cube_m_x in [$t4_min, $t4_max] ")

    t5_min = minimum([t3_min*t4_min, t3_min*t4_max, t3_max*t4_min, t3_max*t4_max])
    t5_max = maximum([t3_min*t4_min, t3_min*t4_max, t3_max*t4_min, t3_max*t4_max])
    println("x1_times_x2 = $x1_times_x2 in [$t5_min, $t5_max]")

    x1_times_x2
end

x0 = [0.75; 0.0; 0.0]
l = [0.5; -1.0; -2.0]
u = [1.0; 1.0; 2.0]
reference_relax(x0, l, u)