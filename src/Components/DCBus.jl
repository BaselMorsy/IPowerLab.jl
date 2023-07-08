using Parameters
@with_kw mutable struct DCBus

    objtyp = 9

    BusID = []
    AreaID = []

    ConnectedLinesIDs = []
    ConnectedLoadsIDs = []
    ConnectedGensIDs  = []
    ConnectedStorageIDs = []


    GeneralSwitch = []

    # single time interval results
    V_magnitude = []

    # multi-time interval results
    V_magnitude_t = []

    # reliability studies results
    V_magnitude_k = []

    # multi-time interval reliability studies results
    V_magnitude_tk = []

    V_max = []
    V_min = []
    V_nom = 1

    P_injected = []

    Pg_net = []
    Pd_net = []

    Coords = []

    measurements = Dict()

end