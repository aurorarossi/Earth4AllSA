module Earth4All

using ModelingToolkit
using WorldDynamics

include("functions.jl")

include("tables.jl")
include("parameters.jl")
include("initialisations.jl")
include("system.jl")

end
