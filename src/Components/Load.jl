using Parameters
@with_kw mutable struct Load

    objtyp = 4

    LoadID = []
    LoadBus_ID = []
    LoadType = [] # To be identefied later on
    
    Pd = []
    Qd = []

    Pd_t = []
    Qd_t = []

    Pd_profile_scaled = []
    Qd_profile_scaled = []

    # single time interval results
    P_shedding = []
    Q_shedding = []

    # multi-time interval results
    P_shedding_t = []
    Q_shedding_t = []

    
    # reliability studies results
    P_shedding_k = []
    Q_shedding_k = []

    
    # multi-time interval reliability studies results
    P_shedding_tk = []
    Q_shedding_tk = []

    Shedding_Cost = 1000
    ActiveControl = false # determines whether the load can be shed or not 

    GeneralSwitch = []
    measurements = []

end