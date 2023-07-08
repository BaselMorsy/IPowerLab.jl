@with_kw mutable struct SubStation

    objtyp = 5

    SubStationBusID = []
    Aux_Buses_IDs = []
    BusbarSections_IDs =[]

    Reconf_AuxLines_IDs = []
    Reconf_CouplerLines_IDs = []

    Branches_IDs = []
    Generators_IDs = []
    Loads_IDs = []
    Storages_IDs = []
    is_split = false
    twins = Dict()

end