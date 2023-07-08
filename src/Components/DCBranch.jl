@with_kw mutable struct DCBranch

    objtyp = 10

    LineID = []
    Fr_bus_ID = []
    To_bus_ID = []

    GeneralSwitch = []

    # single time interval results
    PowerFlow_ij = 0
    PowerFlow_ji = 0
    losses_P = 0

    # multi-time interval results
    PowerFlow_ij_t = 0
    PowerFlow_ji_t = 0
    losses_P_t = 0

    # reliability studies results
    PowerFlow_ij_k = 0
    PowerFlow_ji_k = 0
    losses_P_k = 0

    # multi-time interval reliability studies results
    PowerFlow_ij_tk = 0
    PowerFlow_ji_tk = 0
    losses_P_tk = 0

    r = Inf
    x = Inf
    b = 0

    rating = 50 # pu

    measurements = Dict()

end