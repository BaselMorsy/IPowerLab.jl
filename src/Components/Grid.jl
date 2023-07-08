@with_kw mutable struct PowerGrid
    objtyp = 0
    # Components
    GridID = []
    Areas = []

    Buses = Dict() # AC nodes
    DCBuses = Dict() # DC Nodes

    Branches = Dict() # AC Branches
    DCBranches = Dict() # DC Branches

    Loads = Dict() # AC Loads
    DCLoads = Dict() # DC Loads

    Generators = Dict() # All AC generators including virtual ones
    DCGenerators = Dict() # All DC generators including virtual ones

    Storages = Dict() # AC connected storages
    Substations = Dict() # Switchable substations
    Converters = Dict() # ACDC converters including back to back ACAC converters
    DCLinks = Dict() # DC branches connecting AC nodes
    
    # Necessary helper variables
    Converter_Duplets = Dict()
    Generator_Duplets = Dict()

    # Helper variables that will probably be removed
    Idealized_Branches_Mapping = Dict()
    Idealized_Branches = Dict()
    Idealized_Substations_Mapping = Dict()
    Idealized_Substations = Dict()

    # Internsic Matrices
    Y_bus = []
    Z_bus = []
    b_line = []
    PTDF = []
    LODF = []
    ISF = []

    # Properties
    N_time_steps = 1
    load_profile = []
    generation_profile = []

    N_bus = 0
    N_dc_bus = 0

    N_gen = 0
    N_dc_gen = 0

    N_load = 0
    N_dc_load = 0

    N_branch = 0
    N_dc_branch = 0

    N_substation = 0
    N_switch = 0
    N_storage = 0
    N_converter = 0
    N_dc_links = 0
    
    N_conv_duplets = 0
    N_gen_duplets = 0

    S_base = 100
    load_multiplier = 1
    gen_multiplier = 1
    line_capacity_multiplier = 1
    converter_capacity_multiplier = 1
    Pd_total = []
    Pg_max_total = []

    #Soltuion Queries
    Operating_Cost = []

    BusData_output = DataFrame(Bus=Int64[], V=Float64[], δ=Float64[],
        Pg=Float64[], Qg=Float64[], Pd=Float64[], Qd=Float64[], type=[])

    LineLoading = DataFrame(BranchID=Int64[], FromBus=Int64[], ToBus=Int64[], PL_1=Float64[], PL_2=Float64[],
        PLoss=Float64[], QL_1=Float64[], QL_2=Float64[], QLoss=Float64[], Utilization=Float64[], type=[])

    Converter_flow = DataFrame(ConverterID=Int64[], AC_Bus=Int64[], DC_Bus=Int64[],
        P_ACDC=Float64[], P_DCAC=Float64[], PLoss=Float64[], Utilization=Float64[], type=[])

    Line_Duals = Dict()
    Bus_Duals = Dict()

    z_lines = Dict()
    z_reconf = Dict()
    z_coupler = Dict()

    misc = Dict()
end