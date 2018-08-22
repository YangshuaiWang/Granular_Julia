include("Multi_HZ.jl")
include("Granular_Module.jl")
include("calculate_Pressure.jl")
include("Min_Energy.jl")
include("Deform_Config.jl")
using JuLIP
using Minimise_Energy
using Deform
using Pressure

at = gran.setup_cell(1000, 0.85, 2.0)
min_Energy(at)
at0 = deepcopy(at)
deform_Config(at,"grow", 1.0001)

# include("test.jl")
