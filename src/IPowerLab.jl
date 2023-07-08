module IPowerLab
    using ExportPublic
    using JuMP
    
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