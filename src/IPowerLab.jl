module IPowerLab
    using ExportPublic
    import JuMP, Ipopt
    using JuMP, Ipopt
    import Gurobi
    using Gurobi
    # ENV["GUROBI_HOME"] = "C:\\Program Files\\gurobi903\\win64"
    include("grid_loading.jl")
    include("GridManipulation/grid_functions.jl")
    include("Analysis/PowerMarkets/PowerMarkets.jl")
    include("Analysis/SteadyState/SteadyState.jl")
    include("Parsers/Matpower_Parser.jl")
    include("Parsers/UC_case_parser.jl")
    include("GridVisualization/visualize.jl")
    # include("IntelligentAgents/Jim.jl")

    @exportPublic()
end