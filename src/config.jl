using JuMP, GAMS, SCIP, EAGO, IntervalArithmetic

#=
Creates a model, attaches the appropriate optimizer, 
and sets any relevant parameters
=#

function set_eago_tolerances!(m)
    set_optimizer_attribute(m, "absolute_tolerance", 1E-4)
    set_optimizer_attribute(m, "relative_tolerance", 1E-4)
end

"""
Attaches basic EAGO optimizer (no apriori bilinear relaxation).
"""
function eago()
    m = Model(EAGO.Optimizer)
    set_eago_tolerances!(m)
    set_optimizer_attribute(m, "mul_relax_style", 0)
    return m
end

"""
Attaches EAGO optimizer with subgradient-based bilinear relaxation.
"""
function eago_grad()
    m = Model(EAGO.Optimizer)
    set_eago_tolerances!(m)
    set_optimizer_attribute(m, "mul_relax_style", 1)
    return m
end

"""
Attaches EAGO optimizer with direct enumeration-based bilinear relaxation.
"""
function eago_enum()
    m = Model(EAGO.Optimizer)
    set_eago_tolerances!(m)
    set_optimizer_attribute(m, "mul_relax_style", 3)
    return m
end

"""
Attaches EAGO optimizer with affine arithmetic-based bilinear relaxation.
"""
function eago_affine()
    m = Model(EAGO.Optimizer)
    set_eago_tolerances!(m)
    set_optimizer_attribute(m, "mul_relax_style", 2)
    return m
end

"""
Attaches the baron optimizer from GAMs.
"""
function baron()
    function baron_factory()
        opt = GAMS.Optimizer(GAMS.GAMSWorkspace("C:\\GAMS\\37"))
        opt.gams_options["nlp"] = "BARON"
        opt.gams_options["optca"] = 1E-4
        opt.gams_options["optcr"] = 1E-4
        return opt
    end
    m = Model(baron_factory)
    return m
end

"""
Attaches the scip optimizer.
"""
function scip()
    m = Model(SCIP.Optimizer)
    set_optimizer_attribute(m, "limits/gap", 1E-4)
    set_optimizer_attribute(m, "limits/absgap", 1E-4)
    return m
end
