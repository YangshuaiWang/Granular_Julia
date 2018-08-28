include("Multi_HZ.jl")
include("Granular_Module.jl")
include("calculate_Pressure.jl")
include("Min_Energy.jl")
include("Deform_Config.jl")
using JuLIP
using Minimise_Energy
using Deform
using Pressure

# setup the basic informaiton of the initial config
Grain_Num = 1024
α = 2.0
ϕ_Initial = 0.835

Config_Initial = gran.setup_cell(Grain_Num, ϕ_Initial, α)
while true
    deformType = "grow"
    deformValue = 1.0
    Config_Initial = Deform.deform_Config(Config_Initial, deformType, deformValue)
    V_Initial = get_data(Config_Initial, "Energy")
    if V_Initial < 1.00001 * eps()
        break
    end
end

# setup the growth parameter for fetching the rough jamming configuration
ϕ_Target = 0.845  # this is a volume fraction which is already known bigger than the critical volume fraction
N_Grow = 100      # this is the growth step
ϕ_Increment = (ϕ_Target - ϕ_Initial)/N_Grow
 # the growth history 
Config_Grow_Hist = Array{AbstractAtoms, 1}(N_Grow)
V_Grow_Hist = zeros(N_Grow, 1)
p_Grow_Hist = zeros(N_Grow, 1)

deformType = "phi"
Config_Temp1= deepcopy(Config_Initial)
ϕ_Temp = get_data(Config_Temp1, "Phi")
 # setup the signal for jamming state
Jamming_Signal = 0
for i = 1:N_Grow
    deformValue = ϕ_Temp + ϕ_Increment
    Config_Temp2 = Deform.deform_Config(Config_Temp1, deformType, deformValue)
    Config_Grow_Hist[i] = deepcopy(Config_Temp2)

    V_Grow_Hist[i] = get_data(Config_Temp2, "Energy")
    p_Grow_Hist[i] = trace(get_data(Config_Temp2, "Pressure"))

    if (Jamming_Signal == 0) & (V_Grow_Hist[i] > 5.0e-16) & (p_Grow_Hist[i] > 1e-10)
        Config_RoughJamming = deepcopy(Config_Grow_Hist[i])
        Jamming_Index = i
        Jamming_Signal = 1
    end

    Config_Temp1 = deepcopy(Config_Temp2)
    ϕ_Temp = get_data(Config_Temp1, "Phi")

end
#at = gran.setup_cell(1000, 0.85, 2.0)
#min_Energy(at)
#at0 = deepcopy(at)
#deform_Config(at,"grow", 1.0001)

# include("test.jl")
