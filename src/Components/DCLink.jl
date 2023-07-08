@with_kw mutable struct DCLink

    objtyp = 11

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
    PowerFlow_ij_tk = []
    PowerFlow_ji_tk = []
    losses_P_tk = []

    r = Inf
    x = Inf
    b = 0

    rating = 50 # pu
    Fr_gen_ID = []
    To_gen_ID = []
    measurements = Dict()

end