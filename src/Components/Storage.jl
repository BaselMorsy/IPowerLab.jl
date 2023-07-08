using Parameters
@with_kw mutable struct Storage

    objtyp = 6

    StorageID = []
    StorageBus_ID = []
    StorageType = []

    C0 = []
    C1 = []
    C2 = []

    ramp_up_cost = []
    ramp_down_cost = []

    Pg_max = []
    Pg_min = []
    Qg_max = []
    Qg_min = []

    Pg = []
    Qg = []

    GeneralSwitch = []
    measurements = []

end