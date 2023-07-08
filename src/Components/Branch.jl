using Parameters
@with_kw mutable struct Branch

    objtyp = 2

    LineID = []
    Fr_bus_ID = []
    To_bus_ID = []
    BranchType = 0 # could be 0: for transmission line / transformer, 1: for auxiliary zero impedence line (for substation reconfiguration), or 2: auxiliary zero impedence (busbar coupler)

    GeneralSwitch = []

    # single time interval results
    PowerFlow_ij = 0
    PowerFlow_ji = 0
    losses_P = 0

    ReactFlow_ij = 0
    ReactFlow_ji = 0
    losses_Q = 0

    # multi-time interval results
    PowerFlow_ij_t = 0
    PowerFlow_ji_t = 0
    losses_P_t = 0

    ReactFlow_ij_t = 0
    ReactFlow_ji_t = 0
    losses_Q_t = 0

    # reliability studies results
    PowerFlow_ij_k = 0
    PowerFlow_ji_k = 0
    losses_P_k = 0

    ReactFlow_ij_k = 0
    ReactFlow_ji_k = 0
    losses_Q_k = 0

    # multi-time interval reliability studies results
    PowerFlow_ij_tk = []
    PowerFlow_ji_tk = []
    losses_P_tk = []

    ReactFlow_ij_tk = []
    ReactFlow_ji_tk = []
    losses_Q_tk = []

    r = Inf
    x = Inf
    b = 0

    rating = 50 # pu

    measurements = Dict()

end