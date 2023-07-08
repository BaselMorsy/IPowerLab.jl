using Parameters
@with_kw mutable struct Generator

    objtyp = 3

    GenID = []
    GenBus_ID = []
    GenType = [] # :virtual for p_conv

    C0 = []
    C1 = []
    C2 = []

    ramp_up_cost = []
    ramp_down_cost = []

    start_up_cost = []
    shut_down_cost = []

    min_up_time = []
    min_down_time = []

    Δ_up = []
    Δ_down = []

    Pg_max = []
    Pg_min = []
    Qg_max = []
    Qg_min = []

    # single time interval result
    Pg = []
    Qg = []

    # multi-time interval results
    Pg_t = []
    Qg_t = []

    # reliability studies results
    Pg_k = []
    Qg_k = []

    # multi-time interval reliability studies results
    Pg_tk = []
    Qg_tk = []

    u_gt = []
    α_gt = []
    β_gt = []
    
    GeneralSwitch = []
    Extras = []
    measurements = []

end