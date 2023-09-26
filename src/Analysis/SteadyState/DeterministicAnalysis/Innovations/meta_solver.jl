@with_kw mutable struct MetaSolver
    species = :CCG
    SP_models = Dict()
    MP_models = Dict()
    t_master = []
    t_sub = []
    δ_it = []
    LB_it = []
    UB_it = []
    misc = []
end