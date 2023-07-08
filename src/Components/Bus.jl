using Parameters
@with_kw mutable struct Bus

    objtyp = 1

    BusID = []
    AreaID = []

    ConnectedLinesIDs = []
    ConnectedLoadsIDs = []
    ConnectedGensIDs  = []
    ConnectedStorageIDs = []

    BusType = 0 # Could be 0: non_auxiliary bus, or 1: for auxiliary bus

    GeneralSwitch = []

    # single time interval results
    V_magnitude = []
    δ = []

    # multi-time interval results
    V_magnitude_t = []
    δ_t = []

    # reliability studies results
    V_magnitude_k = []
    δ_k = []

    # multi-time interval reliability studies results
    V_magnitude_tk = []
    δ_tk = []

    V_max = 1.05
    V_min = 0.95
    V_nom = 1

    δ_max = 0.6
    δ_min = -0.6

    P_injected = []

    Pg_net = []
    Pd_net = []

    Qg_net = []
    Qd_net = []

    Coords = []

    measurements = Dict()

end