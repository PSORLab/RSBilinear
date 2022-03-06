using JuMP, GAMS, SCIP, EAGO, IntervalArithmetic

#=
Creates a model, attaches the appropriate optimizer, 
and sets any relevant parameters
=#

function set_eago_tolerances!(m)
    MOI.set(m, MOI.RawOptimizerAttribute("verbosity"), 1)
    MOI.set(m, MOI.RawOptimizerAttribute("cut_min_iterations"), 1)
    MOI.set(m, MOI.RawOptimizerAttribute("cut_max_iterations"), 3)
    MOI.set(m, MOI.RawOptimizerAttribute("output_iterations"), 1)
    MOI.set(m, MOI.RawOptimizerAttribute("iteration_limit"), 100000)
    MOI.set(m, MOI.RawOptimizerAttribute("subgrad_tighten"), true)
    MOI.set(m, MOI.RawOptimizerAttribute("absolute_tolerance"), 1E-6)
    MOI.set(m, MOI.RawOptimizerAttribute("relative_tolerance"), 1E-6)
end

"""
Attaches basic EAGO optimizer (no apriori bilinear relaxation).
"""
function eago()
    opt = EAGO.Optimizer()
    set_eago_tolerances!(opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("mul_relax_style"), 0)
    return opt
end

"""
Attaches EAGO optimizer with subgradient-based bilinear relaxation.
"""
function eago_grad()
    opt = EAGO.Optimizer()
    set_eago_tolerances!(opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("mul_relax_style"), 1)
    return opt
end

"""
Attaches EAGO optimizer with direct enumeration-based bilinear relaxation.
"""
function eago_enum()
    opt = EAGO.Optimizer()
    set_eago_tolerances!(opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("mul_relax_style"), 3)
    return opt
end

"""
Attaches EAGO optimizer with affine arithmetic-based bilinear relaxation.
"""
function eago_affine()
    opt = EAGO.Optimizer()
    set_eago_tolerances!(opt)
    MOI.set(opt, MOI.RawOptimizerAttribute("mul_relax_style"), 2)
    return opt
end

"""
Attaches the baron optimizer from GAMs.
"""
function baron()
    opt = GAMS.Optimizer(GAMS.GAMSWorkspace("C:\\GAMS\\37"))
    opt.gams_options["nlp"] = "BARON"
    opt.gams_options["optca"] = 1E-4
    opt.gams_options["optcr"] = 1E-4
    return opt
end

"""
Attaches the scip optimizer.
"""
function scip()
    opt = SCIP.Optimizer()
    MOI.set(opt, MOI.RawOptimizerAttribute("limits/gap"), 1E-4)
    MOI.set(opt, MOI.RawOptimizerAttribute("limits/absgap"), 1E-4)
    return opt
end
