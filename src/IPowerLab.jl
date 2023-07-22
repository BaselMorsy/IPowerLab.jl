module IPowerLab
    using ExportPublic
    using JuMP, Combinatorics, CSV, JLD2, FileIO, GraphRecipes, Plots, Graphs, GraphPlot
    using JSON, DataStructures, DataFrames, GZip, HTTP
    # using BSON, Statistic, Flux
    
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