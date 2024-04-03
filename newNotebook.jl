### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b286eea2-1258-4cf9-a14a-1eedd94bd387
begin
	using Plots
	using ModelingToolkit
	using DifferentialEquations
	using Statistics
	using PlutoUI
	using IfElse
	using WorldDynamics
	println("Julia environment has been set up")
end

# ╔═╡ b7e9f8cc-f1b1-11ee-19fa-9de56c9de796
md"""
# Making the World Better with Fewer Turnarounds
## A sensitivity analysis of the Earth for All model
"""

# ╔═╡ 0792fd16-759f-4c6d-8973-be41a4510cb8
md"""
	Let's start by setting the Julia environment...
	"""

# ╔═╡ 78c67880-9ba9-410b-be26-27ea66f97d4a
md"""
...and by including the Julia code.
"""

# ╔═╡ fc6c03cf-a831-4967-8446-df9b8d0928d8
begin
_tables = Dict{Symbol,Tuple{Vararg{Float64}}}(
    :ROCWSO => (
        0.06,
        0.02,
        0,
        -0.007,
        -0.01,
    ),
    :IEST => (
        1,
        1,
        0,
    ),
    :PSESTR => (
        0,
        1,
    ),
)
_ranges = Dict{Symbol,Tuple{Float64,Float64}}(
    :IEST => (0, 2),
    :PSESTR => (0, 1),
    :ROCWSO => (0, 2),
)

gettables() = copy(_tables)
getranges() = copy(_ranges)
end;

# ╔═╡ 0c632ed8-0fc0-4f41-a1e8-fb2a40438857
begin
	
#variables: finance - climate - demand - energy - foodland - inventory - labourmarket - other 
_inits = Dict{Symbol,Float64}(
    # climate
    :N2OA => 1.052,
    :CH4A => 2.5,
    :CO2A => 2600,
    :ISCEGA => 12,
    :PWA => 0.4,
    :EHS => 0,
    # demand
    :ETF2022 => 0, # No extra taxes before 2022
    :FGBW => 0.3, # Equal to fraction transfer in 1980
    :POCI => 7081, # It was OCI in 1980
    :WD => 7406.88,
    :PWCIN => 13000,
    :PGCIN => 5400,
    :GNI => 6531.07,
    # energy
    :EEPI2022 => 1,
    :REC => 300,
    :ACSWCF1980 => 10,
    :FEC => 980,
    # finance
    :CCSD => 0.04, # Taken from Vensim table
    :ELTI => 0.02,
    :PI => 0.02,
    :PU => 0.0326951,
    :CBSR => 0.02,
    # foodland
    :BALA => 3000,
    :CRLA => 1450,
    :GRLA => 3300,
    :FOLA => 1100,
    :OGFA => 2600,
    :SQICA => 1,
    :URLA => 215,
    # inventory
    :DELDI => 1,
    :EPP => 28087, # Taken from Vensim table
    :INV => 11234.8, # Taken from Vensim table
    :PRIN => 1, # Taken from Vensim table
    :PRI => 1,
    :RS => 28087,
    :SSWI => 1,
    # labourmarket
    :ECLR => 41,
    :ILPR => 0.8, # Taken from Vensim table
    :LAUS => 3060, # Taken from Vensim table
    :NHW => 2,
    :PURA => 0.05,
    :WARA => 3.6715, # Taken from Vensim table
    :WEOCLR => 1, # Taken from Vensim table
    :WF => 1530,
    :WSO => 0.5,
    # other
    :PGDPP => 6.4 * 0.93,
    # output
    :CUCPIS => 10072.5, # Taken from Vensim table
    :CUCPUS => 909.5, # Taken from Vensim table
    :ETFP => 1,
    :FACNC => 1.05149, # Taken from Vensim table
    :LAUS => 3060, # Taken from Labour and market sector
    :OLY => 26497.1, # Taken from Vensim table
    :WSO => 0.5, # Taken from Labour and market sector
    # population
    :A0020 => 2170,
    :A2040 => 1100,
    :A4060 => 768,
    :A60PL => 382,
    :DEATHS => 30,
    :EGDPP => 6.4,
    :EPA => 0,
    :LE => 67,
    :PA => 62,
    :PASS20 => 100,
    :PASS40 => 64,
    :PASS60 => 38,
    # public
    :RTFPUA => 0, # Taken from Vensim table
    :TFPEE5TA => 1,
    # wellbeing
    :ORP => 0,
    :PAWBI => 0.65,
    :RD => 30,
    :STE => 1.3,
    :SOTR => 0.6,
)

getinitialisations() = copy(_inits)

end;

# ╔═╡ f2d85dc2-d2db-4953-8d49-4051bc5599e8
begin
	_params = Dict{Symbol,Float64}(
    :FOLA1980 => 1100,
    :NHW1980 => 2,
    :WSO1980 => 0.5,
    :WF1980 => 1530,
    :GDPP1980 => 6.4,
    # :INITEGDPP => 6.4,
    :EPA22 => 0,
    # :INITEPA => 0,
    :LE1980 => 67,
    # :INITLE => 67,
    :PA1980 => 62,
    # :INITPA => 62,
    :STE1980 => 1.3,
    :STR1980 => 0.6,
    :ERDN2OKF2022 => 0, # 0.01 in GL,
    :RDN2OKF => 0.01,
    :KN2OKF1980 => 0.11,
    :N2OA1980 => 1.052,
    :MAT => 5,
    :LN2OA => 95,
    :GN2OPP => 5,
    :ERDCH4KC2022 => 0, # 0.01 in GL,
    :RDCH4KC => 0.01,
    :KCH4KC1980 => 0.05,
    :CH4A1980 => 2.5,
    :LCH4A => 7.5,
    :GCH4PP => 5,
    :OBWA2022 => 1.35,
    :SOWLCO2 => 1,
    :LECO2A1980 => 60,
    :CO2A1850 => 2200,
    :TCO2PTCH4 => 2.75,
    :DACCO22100 => 0, # 8 in GL,
    :GCO2PP => 7.9,
    :TCO2ETN2O => 7,
    :TCO2ETCH4 => 23,
    :TCO2ETCO2 => 1,
    :ALGAV => 0.3,
    :ALIS => 0.7,
    :GLSU => 510,
    :ISCEGA1980 => 12,
    :TRSS1980 => 0.01,
    :MRS1980 => 0.0015,
    :WA1980 => 0.4,
    :SVDR => 4,
    :AI1980 => 55,
    :TPM3I => 0.95,
    :HRMI => 333,
    :OWWV => 0.18,
    :WVC1980 => 2,
    :WVF1980 => 0.9,
    :WVWVF => 3,
    :WFEH => 0.0006,
    :EH1980 => 0,
    :PD => 5,
    :TRSA1980 => 0.01,
    :CCCSt => 95,
    :BITRW => 0.2,
    :EETF2022 => 0, # 0.02 in GL,
    :EGTRF2022 => 0, # 0.01 in GL,
    :EPTF2022 => 0, # 0.02 in GL,
    :ETGBW => 0, # 0.2 in GL,
    :FETACPET => 0.5,
    :FETPO => 0.5, # 0.8 in GL,
    :FGDC2022 => 0, # 0.1 in GL,
    :FT1980 => 0.3,
    :GCF => 0.75,
    :GDPOSR => -0.06,
    :GDPP1980 => 6.4,
    :GDDP => 10,
    :GEIC => 0, # 0.02 in GL,
    :GITRO => 0.3,
    :GPP => 200,
    :GSF2022 => 0,
    :INEQ1980 => 0.61,
    :ITRO1980 => 0.4,
    :ITRO2022 => 0.3,
    :OSF1980 => 0.9,
    :MATGF => 0.64,
    :MATWF => 0.39,
    :MWDB => 1,
    :MGDB => 1,
    :STR => 0.03,
    :TAB => 1,
    :TAOC => 1,
    :TAWC => 1,
    :TINT => 5,
    :WCF => 0.9,
    :WDP => 10,
    :WPP => 20,
    :MNFCO2PP => 0.5,
    :FCO2SCCS2022 => 0,
    :GFCO2SCCS => 0.2, # 0.9 in GL,
    :CCCSt => 95,
    :ROCTCO2PT => -0.003,
    :EROCEPA2022 => 0.002, # 0.004 in GL,
    :NIEE => 0.01,
    :GFNE => 0.5, # 1 in GL,
    :FNE2022 => 0.03,
    :FNE1980 => 0,
    :EUEPRUNEFF => 3,
    :ECRUNEFF => 10,
    :GREF => 0.5, # 1 in GL,
    :REFF2022 => 0.23,
    :REFF1980 => 0.065,
    :RCUT => 3,
    :RECT => 3,
    :LREC => 40,
    :SWC1980 => 10,
    :CRDSWC => 0.2,
    :CAPEXRE1980 => 7,
    :CAPEXFED => 0.7,
    :OPEXRED => 0.001,
    :OPEXFED => 0.02,
    :CNED => 0.033,
    :FREH => 0,
    :KWEPKGH2 => 40,
    :TPTH2 => 10,
    :BEM => 0,
    :EFPP => 0.345,
    :TWHPEJCE => 278,
    :MTPEJCE => 24,
    :EKHPY => 8,
    :FECCT => 3,
    :NLFEC => 40,
    :sFCUTLOFC => 0.5,
    :NCUT => 8,
    :TCE => 0.03,
    :AFMCM => 1.35,
    :TCFFFNEU => 240,
    :TC => 0.02,
    :FSRT => 1,
    :GRCR => 0,
    :IEFT => 10,
    :INSR => 0.7,
    :IPTCB => 1,
    :IT => 0.02,
    :NBBM => 0.005,
    :NBOM => 0.015,
    :NSR => 0.02,
    :SRAT => 1,
    :UNSR => -1.5,
    :UPTCB => 1,
    :UT => 0.05,
    :AFGDP => 0.05,
    :CBECLE => -0.03,
    :CO2ARA => 1,
    :CO2C2022 => 420,
    :CO2CEACY => 0.3,
    :CO2RHFC => 65,
    :CRDRA => 0.05,
    :CTF => 500,
    :CYRA => 5,
    :DRC => 0.05,
    :ECRA22 => 400,
    :EGB22 => 5,
    :EROCFSP => 0,
    :FCG => 0.1,
    :FFLREOGRRM => -5,
    :FU80 => 61,
    :FF80 => 88450,
    :FUELER => 0.02,
    :FUESQ => -0.001,
    :GCWR => 0.05, # 0.2 in GL,
    :GFNRM => 0.1, # 0.5 in GL,
    :GFRA => 0.1, # 0.5 in GL,
    :KCKRM => 24,
    :LER80 => 0.004,
    :MFAM => 2,
    :OGRR80 => 0.004,
    :OW2022 => 1.35,
    :OWEACY => -0.3,
    :ROCFP => 0.01,
    :ROCFSP => 0.002,
    :SFU => 20,
    :SSP2LMA => 1,
    :TCTB => 0,
    :TFFLR => 0.2,
    :UDT => 10,
    :ULP => 0.05,
    :DAT => 1.2,
    :DDI1980 => 1,
    :DIC => 0.4,
    :DRI => 1,
    :ICPT => 0.25,
    :MRIWI => 1.07,
    :OO => 28087,
    :PH => 0,
    :PPU => 1,
    :SAT => 1,
    :INVEODDI => -0.6,
    :INVEOIN => -0.26,
    :INVEOSWI => -0.6,
    :SRI => 1,
    :SWI1980 => 1,
    :TAS => 0.24,
    :AUR => 0.05,
    :FIC => 1,
    :GDPP1980 => 6.4, # Should be the same in Population
    :GDPPEROCCLRM => -0.1,
    :GENLPR => 0,
    :NLPR80 => 0.85,
    :PFTJ => 1,
    :PFTJ80 => 1,
    :PRUN => 1, # Taken from Vensim table
    :PUELPR => 0.05,
    :ROCECLR80 => 0.02,
    :RWER => 0.015,
    :TAHW => 5,
    :TENHW => -0.03,
    :TELLM => 5,
    :TYLD => 2.3,
    :WSOECLR => 1.05,
    :WSOELPR => 0.2,
    :TEGR => 4,
    :INELOK => -0.5,
    :NK => 0.3,
    :CAPPIS1980 => 59250,
    :CAPPUS1980 => 5350,
    :CC1980 => 1,
    :CTPIS => 1.5,
    :CTPUS => 1.5,
    :EMCUC => 1.7,
    :FCI => 0,
    :LAUS1980 => 3060, # Taken from Vensim table
    :LCPIS1980 => 15,
    :OW2022 => 1.35, # Taken from Climate sector
    :OWECCM => 0.2,
    :OWELCM => -0.1,
    :USPIS2022 => 0, # 0.01 in GL,
    :USPUS2022 => 0, # 0.01 in GL,
    :CBCEFRA => -0.8,
    :CU1980 => 0.8,
    :ED1980 => 1,
    :EDEFRA => 5,
    :EDELCM => 0.5,
    :FRA1980 => 0.9,
    :FRACAM => 0.65,
    :GDPP1980 => 6.4, # It should be the same as in the Labour and market sector
    :GDPPEFRACA => -0.2,
    :IPT => 1,
    :JOBS1980 => 1600,
    :KAPPA => 0.3,
    :LAMBDA => 0.7, # Calculated as 1-KAPPA
    :MA1980 => 0.25,
    :OG1980 => 0.06,
    :OO1980 => 28087, # Taken from Vensim table
    :PCORPIS => 2.3,
    :PCORPUS => 2.3,
    :TOED => 1,
    :WSOEFRA => -2.5,
    :CMFR => 0.01,
    :DNC80 => 4.3,
    :DNCA => 0,
    :DNCG => 0.14,
    :DNCM => 1.2,
    :EIP => 30,
    :EIPF => 0,
    :FP => 20,
    :FW => 0.5,
    :FADFS => 0.8,
    :GEFR => 0, # 0.2 in GL,
    :GEPA => 0,
    :LEA => 0.001,
    :LEEPA => 0.75,
    :LEG => 0.15,
    :LEMAX => 85,
    :MFM => 1.6,
    :MLEM => 1.1,
    :ORDER => 10,
    :OW2022 => 1.35,
    :OWELE => -0.02,
    :SSP2FA2022F => 1,
    :TAHI => 10,
    :CTA2022 => 9145,
    :CTPIS => 1.5, # Taken from Output sector
    :EDROTA2022 => 0.003,
    :DROTA1980 => 0.01,
    :FUATA => 0.3,
    :GDPTL => 15,
    :IIEEROTA => -0.1,
    :IPR1980 => 1.2,
    :IPRVPSS => 1,
    :IPT => 1, # Taken from Output sector
    :MIROTA2022 => 0, # 0.005 in GL,
    :OWETFP => -0.1,
    :SC1980 => 0.3,
    :SCROTA => 0.5,
    :XETAC2022 => 0,
    :XETAC2100 => 0,
    :OW2022 => 1.35,
    :AI => 0.6,
    :AP => 0.02,
    :AWBPD => 9,
    :DRDI => 0.5,
    :DRPS => 0.7,
    :EIP => 30,
    :EIPF => 0,
    :GWEAWBGWF => -0.58,
    :IEAWBIF => -0.6,
    :MWBGW => 0.2,
    :NRD => 30,
    :PESTF => -15,
    :PAEAWBF => 0.5,
    :PREAWBF => 6,
    :SPS => 0.3,
    :STEERDF => 1,
    :STRERDF => -1,
    :TCRD => 10,
    :TDI => 15,
    :TEST => 10,
    :TI => 0.5,
    :THPA => 0.8,
    :TPR => 0.02,
    :TPS => 3,
    :TW => 1,
)

getparameters() = copy(_params)

end;

# ╔═╡ a30de87a-fae9-4270-b1da-3a80de97b235
begin


	function add_equation!(eqs, equation)
   append!(eqs, [equation])
end

function delay_n!(eqs, x, rt, lv, delay, order)
   append!(eqs, [rt[1] ~ lv[1] / (delay / order)])
   append!(eqs, [D(lv[1]) ~ x - rt[1]])
   for d in 2:order
      append!(eqs, [rt[d] ~ lv[d] / (delay / order)])
      append!(eqs, [D(lv[d]) ~ rt[d-1] - rt[d]])
   end
end

ramp(x, slope, startx, endx) = IfElse.ifelse(x > startx, IfElse.ifelse(x < endx, slope * (x - startx), slope * (endx - startx)), 0)

function smooth!(eqs, x, input, delay_time)
   append!(eqs, [D(x) ~ (input - x) / delay_time])
end

clip(returnifgte, returniflt, inputvalue, threshold) = IfElse.ifelse(inputvalue ≥ threshold, returnifgte, returniflt)

step(inputvalue, returnifgte, threshold) = clip(returnifgte, zero(returnifgte), inputvalue, threshold)


pulse(inputvalue, start, width) = IfElse.ifelse(inputvalue >= start, 1, 0) * IfElse.ifelse(inputvalue < (start + width), 1, 0)

interpolate(x, x₁, xₙ, y₁, yₙ) = y₁ + (x - x₁) * ((yₙ - y₁) / (xₙ - x₁))

function interpolate(x, yvalues::Tuple{Vararg{Float64}}, xrange::Tuple{Float64,Float64})
   interpolate(x, collect(yvalues), collect(LinRange(xrange[1], xrange[2], length(yvalues))))
end

function interpolate(x, yvalues::Vector{Float64}, xvalues::Vector{Float64})
   # y gets the min y value if less than the min x range
   #   or the max y value if greater than the max x range
   y = (x < xvalues[1]) * yvalues[1] + (x ≥ xvalues[end]) * yvalues[end]

   # in case x is inside the range, y gets the interpolated value
   for i ∈ 1:length(yvalues)-1
      y += (x ≥ xvalues[i]) * (x < xvalues[i+1]) * interpolate(x, xvalues[i], xvalues[i+1], yvalues[i], yvalues[i+1])
   end

   return y
end

function withlookup(x, pairs::Vector{Tuple{Float64,Float64}})
   interpolate(x, map(t -> t[end], pairs), map(t -> t[1], pairs))
end
	
	@variables t
	D = Differential(t)
	
function earth4all(; name, params=_params, inits=_inits, tables=_tables, ranges=_ranges)
    ORDER = Int64(params[:ORDER])

    eqs = []

    # climate
    CCCSt = params[:CCCSt]
    @parameters ERDN2OKF2022 = params[:ERDN2OKF2022] [description = "Extra rate of decline in N2O per kg fertilizer from 2022"]
    RDN2OKF = params[:RDN2OKF]
    KN2OKF1980 = params[:KN2OKF1980]
    N2OA1980 = params[:N2OA1980]
    MAT = params[:MAT]
    LN2OA = params[:LN2OA]
    GN2OPP = params[:GN2OPP]
    @parameters ERDCH4KC2022 = params[:ERDCH4KC2022] [description = "Extra rate of decline in CH4 per kg crop after 2022 1/y"]
    RDCH4KC = params[:RDCH4KC]
    KCH4KC1980 = params[:KCH4KC1980]
    CH4A1980 = params[:CH4A1980]
    LCH4A = params[:LCH4A]
    GCH4PP = params[:GCH4PP]
    OBWA2022 = params[:OBWA2022]
    SOWLCO2 = params[:SOWLCO2]
    LECO2A1980 = params[:LECO2A1980]
    CO2A1850 = params[:CO2A1850]
    TCO2PTCH4 = params[:TCO2PTCH4]
    @parameters DACCO22100 = params[:DACCO22100] [description = "Direct air capture of CO2 in 2100 GtCO2/y"]
    GCO2PP = params[:GCO2PP]
    TCO2ETN2O = params[:TCO2ETN2O]
    TCO2ETCH4 = params[:TCO2ETCH4]
    TCO2ETCO2 = params[:TCO2ETCO2]
    ALGAV = params[:ALGAV]
    ALIS = params[:ALIS]
    GLSU = params[:GLSU]
    ISCEGA1980 = params[:ISCEGA1980]
    TRSS1980 = params[:TRSS1980]
    MRS1980 = params[:MRS1980]
    WA1980 = params[:WA1980]
    SVDR = params[:SVDR]
    AI1980 = params[:AI1980]
    TPM3I = params[:TPM3I]
    HRMI = params[:HRMI]
    OWWV = params[:OWWV]
    WVC1980 = params[:WVC1980]
    WVF1980 = params[:WVF1980]
    WVWVF = params[:WVWVF]
    WFEH = params[:WFEH]
    EH1980 = params[:EH1980]
    TRSA1980 = params[:TRSA1980]
    PD = params[:PD]

    @variables AL(t) [description = "ALbedo (1)"]
    @variables AL1980(t) [description = "ALbedo in 1980 (1)"]
    @variables CAC(t) [description = "Cost of air capture GDollar/y"]
    @variables CH4A(t) = inits[:CH4A] [description = "CH4 in Atmosphere GtCH4"]
    @variables CH4BD(t) [description = "CH4 BreakDown GtCH4/y"]
    @variables CH4CA(t) [description = "CH4 concentration in atm ppm"]
    @variables CH4E(t) [description = "CH4 emissions GtCH4/y"]
    @variables CH4FPP(t) [description = "CH4 forcing per ppm W/m2/ppm"]
    @variables CO2A(t) = inits[:CO2A] [description = "CO2 in Atmosphere GtCO2"]
    @variables CO2AB(t) [description = "CO2 absorption GtCO2/y"]
    @variables CO2CA(t) [description = "CO2 concentration in atm ppm"]
    @variables CO2E(t) [description = "CO2 emissions GtCO2/y"]
    @variables CO2FCH4(t) [description = "CO2 from CH4 GtCO2/y"]
    @variables CO2FPP(t) [description = "CO2 forcing per ppm W/m2/ppm"]
    @variables CO2GDP(t) [description = "CO2 per GDP (kgCO2/Dollar)"]
    @variables DACCO2(t) [description = "Direct Air Capture of CO2 GtCO2/y"]
    @variables ECIM(t) [description = "Extra cooling from ice melt ZJ/y"]
    @variables EHS(t) = inits[:EHS] [description = "Extra heat in surface ZJ"]
    @variables EWFF(t) [description = "Extra Warming from forcing ZJ/y"]
    @variables FCH4(t) [description = "Forcing from CH4 W/m2"]
    @variables FCO2(t) [description = "Forcing from CO2 W/m2"]
    @variables FN2O(t) [description = "Forcing from N2O W/m2"]
    @variables FOG(t) [description = "Forcing from other gases W/m2"]
    @variables GHGE(t) [description = "GHG emissions GtCO2e/y"]
    @variables HDO(t) [description = "Heat to deep ocean ZJ/y"]
    @variables HTS(t) [description = "Heat to space ZJ/y"]
    @variables ISC(t) [description = "Ice and snow cover Mha"]
    @variables ISCEGA(t) = inits[:ISCEGA] [description = "Ice and snow cover excl G&A Mkm2"]
    @variables KCH4EKC(t) [description = "kg CH4 emission per kg crop"]
    @variables KN2OEKF(t) [description = "kg N2O emission per kg fertiliser"]
    @variables LECO2A(t) [description = "Life of extra CO2 in atm y"]
    @variables MEL(t) [description = "Melting Mha/y"]
    @variables MMCH4E(t) [description = "Man-made CH4 emissions GtCH4/y"]
    @variables MMF(t) [description = "Man-made Forcing W/m2"]
    @variables MMN2OE(t) [description = "Man-made N2O emissions GtN2O/y"]
    @variables MRDI(t) [description = "Melting rate deep ice 1/y"]
    @variables MRS(t) [description = "Melting rate surface 1/y"]
    @variables N2OA(t) = inits[:N2OA] [description = "N2O in Atmosphere GtN2O"]
    @variables N2OBD(t) [description = "N2O BreakDown GtN2O/y"]
    @variables N2OCA(t) [description = "N2O concentration in atm ppm"]
    @variables N2OE(t) [description = "N2O emissions GtN2O/y"]
    @variables N2OFPP(t) [description = "N2O forcing per ppm W/m2/ppm"]
    @variables NCH4E(t) [description = "Natural CH4 emissions GtCH4/y"]
    @variables NN2OE(t) [description = "Natural N2O emissions GtN2O/y"]
    @variables OW(t) [description = "OBserved WArming deg C"]
    @variables OWLCO2(t) [description = "OWeoLoCO2"]
    @variables PWA(t) = inits[:PWA] [description = "Perceived WArming deg C"]
    @variables REHE(t) [description = "Risk of extreme heat event (1)"]
    @variables TMMF(t) [description = "Total man-made forcing W/m2"]
    @variables TRHGA(t) [description = "Transfer rate for heat going to abyss 1/y"]
    @variables TRHGS(t) [description = "Transfer rate for heat going to space 1/y"]
    @variables WVC(t) [description = "Water Vapor Concentration g/kg"]
    @variables WVF(t) [description = "Water Vapour Feedback W/m2"]

    # demand
    BITRW = params[:BITRW]
    @parameters EETF2022 = params[:EETF2022] [description = "Extra Empowerment Tax From 2022 (share of NI)"]
    @parameters EGTRF2022 = params[:EGTRF2022] [description = "Extra General Tax Rate From 2022"]
    @parameters EPTF2022 = params[:EPTF2022] [description = "Extra Pension Tax From 2022 (share of NI)"]
    @parameters ETGBW = params[:ETGBW] [description = "Extra Transfer of Govmnt Budget to Workers"]
    FETACPET = params[:FETACPET]
    @parameters FETPO = params[:FETPO] [description = "Fraction of Extra Taxes Paid by Owners"]
    @parameters FGDC2022 = params[:FGDC2022] [description = "Fraction of Govmnt Debt Cancelled in 2022 1/y"]
    FT1980 = params[:FT1980]
    GCF = params[:GCF]
    GDDP = params[:GDDP]
    GDPOSR = params[:GDPOSR]
    GDPP1980 = params[:GDPP1980]
    @parameters GEIC = params[:GEIC] [description = "Goal for Extra Income from Commons (share of NI)"]
    GITRO = params[:GITRO]
    GPP = params[:GPP]
    GSF2022 = params[:GSF2022]
    INEQ1980 = params[:INEQ1980]
    ITRO1980 = params[:ITRO1980]
    ITRO2022 = params[:ITRO2022]
    MATGF = params[:MATGF]
    MATWF = params[:MATWF]
    MGDB = params[:MGDB]
    MWDB = params[:MWDB]
    OSF1980 = params[:OSF1980]
    STR = params[:STR]
    TAB = params[:TAB]
    TAOC = params[:TAOC]
    TAWC = params[:TAWC]
    TINT = params[:TINT]
    WCF = params[:WCF]
    WDP = params[:WDP]
    WPP = params[:WPP]

    @variables BCIL(t) [description = "Bank Cash Inflow from Lending GDollar/y"]
    @variables BCISNI(t) [description = "Bank Cash Inflow as Share of NI (1)"]
    @variables BITRO(t) [description = "Basic Income Tax Rate Owners (1)"]
    @variables CANCD(t) [description = "CANCellation of Debt GDollar/y"]
    @variables CD(t) [description = "Consumption Demand GDollar/y"]
    @variables CFGB(t) [description = "Cash Flow from Govmnt to Banks GDollar/y"]
    @variables CFWB(t) [description = "Cash Flow from Workers to Banks GDollar/y"]
    @variables CONTR(t) [description = "CONTRol: (C+G+S)/NI = 1"]
    @variables CPP(t) [description = "Consumption Per Person GDollar/y"]
    @variables CSGDP(t) [description = "Consumption Share of GDP (1)"]
    @variables EGTF2022(t) [description = "Extra General Tax From 2022 Gdollar/y"]
    @variables ETF2022(t) = inits[:ETF2022] [description = "Extra Taxes From 2022 GDollar/y"]
    @variables ETTAF2022(t) [description = "Extra Taxes for TAs From 2022 GDollar/y"]
    @variables FGBW(t) = inits[:FGBW] [description = "Fraction of Govmnt Budget to Workers (1)"]
    @variables GCIN(t) [description = "Govmnt Cash INflow GDollar/y"]
    @variables GD(t) = 28087 * params[:MATGF] [description = "Govmnt Debt Gdollar"]
    @variables GDB(t) [description = "Govmnt Debt Burden y"]
    @variables GETF2022(t) [description = "Goal for Extra Taxes From 2022 GDollar/y"]
    @variables GFGBW(t) [description = "Goal for Fraction of Govmnt Budget to Workers (1)"]
    @variables GFSNI(t) [description = "Govmnt Finance as Share of NI (1)"]
    @variables GGI(t) [description = "Govmnt Gross Income GDollar/y"]
    @variables GGIS(t) [description = "Govmnt Gross Income (as Share of NI)"]
    @variables GIC(t) [description = "Govmnt Interest Cost GDollar/y"]
    @variables GIPC(t) [description = "Govmnt Investment in Public Capacity GDollar/y"]
    @variables GND(t) [description = "Govmnt New Debt GDollar/y"]
    @variables GNI(t) = inits[:GNI] [description = "Govmnt net income GDollar/y"]
    @variables GNISNI(t) [description = "Govmnt Net Income as Share of NI (1)"]
    @variables GP(t) [description = "Govmnt Payback GDollar/y"]
    @variables GPU(t) [description = "Govmnt PUrchases GDollar/y"]
    @variables GS(t) [description = "Govmnt Spending GDollar/y"]
    @variables GSGDP(t) [description = "Govmnt Share of GDP (1)"]
    @variables IC2022(t) [description = "Income from Commons from 2022 GDollar/y"]
    @variables INEQ(t) [description = "INEQuality (1)"]
    @variables INEQI(t) [description = "INEQuality Index (1980=1)"]
    @variables ITO(t) [description = "Income Tax Owners (1)"]
    @variables ITW(t) [description = "Income Tax Workers (1)"]
    @variables MGD(t) [description = "Max Govmnt Debt GDollar"]
    @variables MWD(t) [description = "Max Workers Debt GDollar"]
    @variables OC(t) [description = "Owner Consumption GDollar/y"]
    @variables OCF(t) [description = "Owner Consumptin Fraction (1)"]
    @variables OCIN(t) [description = "Owner Cash INflow GDollar/y"]
    @variables OI(t) [description = "Owner Income GDollar/y"]
    @variables OOIAT(t) [description = "Owner Operating Income After Tax GDollar/y"]
    @variables OS(t) [description = "Owner Savings GDollar/y"]
    @variables OSF(t) [description = "Owner Savings Fraction (1)"]
    @variables OT(t) [description = "Owner Taxes GDollar/y"]
    @variables OTR(t) [description = "Owner Tax Rate (1)"]
    @variables PGCIN(t) = inits[:PGCIN] [description = "Permanent govmnt cash inflow GDollar/y"]
    @variables POCI(t) = inits[:POCI] [description = "Permanent Owner Cash Inflow GDollar/y"]
    @variables PWCIN(t) = inits[:PWCIN] [description = "Permanent Worker Cash INflow GDollar/y"]
    @variables SSGDP(t) [description = "Savings Share of GDP (1)"]
    @variables ST(t) [description = "Sales Tax GDollar/y"]
    @variables STO(t) [description = "Sales Tax Owners GDollar/y"]
    @variables STW(t) [description = "Sales Tax Workers GDollar/y"]
    @variables TP(t) [description = "Transfer Payments GDollar/y"]
    @variables TPP(t) [description = "Total Purchasing Power GDollar/y"]
    @variables TS(t) [description = "Total Savings GDollar/y"]
    @variables WCD(t) [description = "Worker consumption demand GDollar/y"]
    @variables WCIN(t) [description = "Worker Cash INflow GDollar/y"]
    @variables WD(t) = 18992 * params[:MATWF] [description = "Workers Debt GDollar"]
    @variables WDB(t) [description = "Worker Debt Burden y"]
    @variables WDI(t) [description = "Worker Disposable Income kDollar/p/y"]
    @variables WFCSI(t) [description = "Worker Finance Cost as Share of Income (1)"]
    @variables WI(t) [description = "Worker Income GDollar/y"]
    @variables WIAT(t) [description = "Worker Income After Tax GDollar/y"]
    @variables WIC(t) [description = "Worker Interest Cost GDollar/y"]
    @variables WND(t) [description = "Workers New Debt GDollar/y"]
    @variables WP(t) [description = "Workers Payback GDollar/y"]
    @variables WS(t) [description = "Worker Savings GDollar/y"]
    @variables WT(t) [description = "Worker Taxes GDollar/y"]
    @variables WTR(t) [description = "Worker Tax Rate (1)"]

    # energy
    MNFCO2PP = params[:MNFCO2PP]
    FCO2SCCS2022 = params[:FCO2SCCS2022]
    @parameters GFCO2SCCS = params[:GFCO2SCCS] [description = "Goal for fraction of CO2-sources with CCS"]
    CCCSt = params[:CCCSt]
    ROCTCO2PT = params[:ROCTCO2PT]
    @parameters EROCEPA2022 = params[:EROCEPA2022] [description = "Extra ROC in energy productivity after 2022 1/y"]
    NIEE = params[:NIEE]
    @parameters GFNE = params[:GFNE] [description = "Goal for fraction new electrification"]
    FNE2022 = params[:FNE2022]
    FNE1980 = params[:FNE1980]
    EUEPRUNEFF = params[:EUEPRUNEFF]
    ECRUNEFF = params[:ECRUNEFF]
    @parameters GREF = params[:GREF] [description = "Goal for renewable el fraction"]
    REFF2022 = params[:REFF2022]
    REFF1980 = params[:REFF1980]
    RCUT = params[:RCUT]
    RECT = params[:RECT]
    LREC = params[:LREC]
    SWC1980 = params[:SWC1980]
    CRDSWC = params[:CRDSWC]
    CAPEXRE1980 = params[:CAPEXRE1980]
    CAPEXFED = params[:CAPEXFED]
    OPEXRED = params[:OPEXRED]
    OPEXFED = params[:OPEXFED]
    CNED = params[:CNED]
    FREH = params[:FREH]
    KWEPKGH2 = params[:KWEPKGH2]
    TPTH2 = params[:TPTH2]
    BEM = params[:BEM]
    EFPP = params[:EFPP]
    TWHPEJCE = params[:TWHPEJCE]
    MTPEJCE = params[:MTPEJCE]
    EKHPY = params[:EKHPY]
    FECCT = params[:FECCT]
    NLFEC = params[:NLFEC]
    sFCUTLOFC = params[:sFCUTLOFC]
    NCUT = params[:NCUT]
    TCE = params[:TCE]
    AFMCM = params[:AFMCM]
    TCFFFNEU = params[:TCFFFNEU]
    TC = params[:TC]

    @variables CNEL(t) [description = "Cost of Nuclear ELectricity Gdollar/y"]
    @variables CO2NFIP(t) [description = "CO2 from non-fossil industrial processes GtCO2/y"]
    @variables FCO2SCCS(t) [description = "Fraction of CO2-sources with CCS (1)"]
    @variables NFCO2PP(t) [description = "Non-fossil CO2 per person tCO2/p/y"]
    @variables CO2EI(t) [description = "CO2 from energy and industry GtCO2/y"]
    @variables CO2EP(t) [description = "CO2 from energy production GtCO2/y"]
    @variables CO2EMPP(t) [description = "CO2 EMissions per person tCO2/y"]
    @variables CCCSG(t) [description = "Cost of CCS GDollar/y"]
    @variables ICCSC(t) [description = "Installed CCS capacity GtCO2/y"]
    @variables TCO2PT(t) [description = "tCO2 per toe"]
    @variables EEPI2022(t) = inits[:EEPI2022] [description = "Extra energy productivity index 2022=1"]
    @variables IEEPI(t) [description = "Increase in extra energy productivity index 1/y"]
    @variables TPPUEBEE(t) [description = "Traditional per person use of electricity before EE MWh/p/y"]
    @variables TPPUFFNEUBEE(t) [description = "Traditional per person use of fossil fuels for non-el-use before EE toe/p/y"]
    @variables DEBNE(t) [description = "Demand for electricity before NE TWh/y"]
    @variables DFFNEUBNE(t) [description = "Demand for fossil fuel for non-el use before NE Mtoe/y"]
    @variables FNE(t) [description = "Fraction new electrification (1)"]
    @variables ERDNEFFFNE(t) [description = "Extra reduction in demand for non-el fossil fuel from NE Mtoe/y"]
    @variables CNE(t) [description = "Cost of new electrification GDollar/y"]
    @variables EIDEFNE(t) [description = "Extra increase in demand for electricity from NE TWh/y"]
    @variables DFFFNEU(t) [description = "Demand for fossil fuel for non-el use Mtoe/y"]
    @variables UFF(t) [description = "Use of fossil fuels Mtoe/y"]
    @variables DE(t) [description = "Demand for electricity TWh/y"]
    @variables DRES(t) [description = "Desired renewable electricity share (1)"]
    @variables DSRE(t) [description = "Desired supply of renewable electricity TWh/y"]
    @variables DREC(t) [description = "Desired renewable el capacity GW"]
    @variables DRECC(t) [description = "Desired renewable el capacity change GW"]
    @variables REC(t) = inits[:REC] [description = "Renewable electricity capacity GW"]
    @variables AREC(t) [description = "Addition of renewable el capacity GW/y"]
    @variables DIREC(t) [description = "Discard of renewable el capacity GW/y"]
    @variables ASWC(t) [description = "Addition of sun and wind capacity GW/y"]
    @variables ACSWCF1980(t) = inits[:ACSWCF1980] [description = "ACcumulated sun and wind capacity from 1980 GW"]
    @variables NDSWC(t) [description = "Number of doublings in sun and wind capacity (1)"]
    @variables CISWC(t) [description = "Cost index for sun and wind capacity (1)"]
    @variables CAPEXRED(t) [description = "CAPEX renewable el dollar/W"]
    @variables CAPEXREG(t) [description = "CAPEX renewable el Gdollar/y"]
    @variables OPEXREG(t) [description = "OPEX renewable el Gdollar/y"]
    @variables CRE(t) [description = "Cost of renewable electricity GDollar/y"]
    @variables CAPEXFEG(t) [description = "CAPEX fossil el GDollar/y"]
    @variables OPEXFEG(t) [description = "OPEX fossil el GDollar/y"]
    @variables CFE(t) [description = "Cost of fossil electricity GDollar/y"]
    @variables CEL(t) [description = "Cost of electricity GDollar/y"]
    @variables REP(t) [description = "Renewable electricity production TWh/y"]
    @variables GHMH2(t) [description = "Green hydrogen MtH2/y"]
    @variables GHMt(t) [description = "Green hydrogen Mtoe/y"]
    @variables RHP(t) [description = "Renewable heat production Mtoe/y"]
    @variables TWEPEJEE(t) [description = "TWh-el per EJ - engineering equivalent"]
    @variables IIASAREP(t) [description = "IIASA Renewable energy production EJ/yr"]
    @variables FTWEPMt(t) [description = "4 TWh-el per Mtoe"]
    @variables IIASAFEP(t) [description = "IIASA Fossil energy production EJ/yr"]
    @variables LCEP(t) [description = "Low-carbon el production TWh/y"]
    @variables DFE(t) [description = "Demand for fossil electricity TWh/y"]
    @variables DFEC(t) [description = "Desired fossil el capacity GW"]
    @variables DFECC(t) [description = "Desired fossil el capacity change GW/y"]
    @variables AFEC(t) [description = "Addition of fossil el capacity GW/y"]
    @variables FEC(t) = inits[:FEC] [description = "Fossil electricity capacity GW"]
    @variables LFEC(t) [description = "Life of fossil el capacity y"]
    @variables DIFEC(t) [description = "DIscard of Fossil El Capacity GW/y"]
    @variables FCUT(t) [description = "Fossil capacity up-time kh/y"]
    @variables FCUTLOFC(t) [description = "FCUTeoLOFC (1)"]
    @variables FEP(t) [description = "Fossil Electricity Production TWh/y"]
    @variables NC(t) [description = "Nuclear capacity GW"]
    @variables NEP(t) [description = "Nuclear electricity production TWh/y"]
    @variables EP(t) [description = "Electricity production TWh/y"]
    @variables ELB(t) [description = "ELectricity balance (1)"]
    @variables FFPNE(t) [description = "Fraction fossil plus nuclear electricity (1)"]
    @variables EU(t) [description = "Energy use Mtoe/y"]
    @variables EUPP(t) [description = "Energy use per person toe/p/y"]
    @variables FFE(t) [description = "Fossil fuels for electricity Mtoe/y"]
    @variables TCEG(t) [description = "Traditional cost of electricity GDollar/y"]
    @variables TCFFFNEUG(t) [description = "Traditional cost of fossil fuel for non-el use Gdollar/y"]
    @variables CFFFNEU(t) [description = "Cost of Fossil Fuel For Non-El Use Gdollar/y"]
    @variables CG(t) [description = "Cost of grid GDollar/y"]
    @variables TGC(t) [description = "Traditional grid cost GDollar/y"]
    @variables TCEN(t) [description = "Traditional Cost of ENergy Gdollar/y"]
    @variables TCENSGDP(t) [description = "Traditional cost of energy as share of GDP (1)"]
    @variables CE(t) [description = "Cost of energy GDollar/y"]
    @variables RECTEC(t) [description = "Ratio of Energy cost to Trad Energy cost (1)"]
    @variables CESGDP(t) [description = "Cost of energy as share of GDP (1)"]
    @variables ECETSGDP(t) [description = "Extra cost of Energy Turnaround as share of GDP (1)"]

    # finance
    FSRT = params[:FSRT]
    GRCR = params[:GRCR]
    IEFT = params[:IEFT]
    INSR = params[:INSR]
    IPTCB = params[:IPTCB]
    IT = params[:IT]
    NBOM = params[:NBOM]
    NBBM = params[:NBBM]
    NSR = params[:NSR]
    SRAT = params[:SRAT]
    UPTCB = params[:UPTCB]
    UNSR = params[:UNSR]
    UT = params[:UT]

    @variables CBC(t) [description = "Corporate Borrowing Cost 1/y"]
    @variables CBSR(t) = inits[:CBSR] [description = "Central Bank Signal Rate 1/y"]
    @variables CCSD(t) = inits[:CCSD] [description = "Cost of Capital for Secured Debt 1/y"]
    @variables CSR(t) [description = "Change in Signal Rate 1/yy"]
    @variables GBC(t) [description = "Govmnt Borrowing Cost 1/y"]
    @variables CBC1980(t) [description = "Corporate Borrowing Cost in 1980 1/y"]
    @variables ELTI(t) = inits[:ELTI] [description = "Expected Long Term Inflation 1/y"]
    @variables ISR(t) [description = "Indicated Signal Rate 1/y"]
    @variables NCCR(t) [description = "Normal Corporate Credit Risk 1/y"]
    @variables PI(t) = inits[:PI] [description = "Perceived Inflation CB 1/y"]
    @variables PU(t) = inits[:PU] [description = "Perceived Unemployment CB (1)"]
    @variables TGIR(t) [description = "10-yr Govmnt Interest Rate 1/y"]
    @variables TIR(t) [description = "3m Interest Rate 1/y"]
    @variables WBC(t) [description = "Worker Borrowing Cost 1/y"]

    # foodland
    AFGDP = params[:AFGDP]
    CBECLE = params[:CBECLE]
    CO2ARA = params[:CO2ARA]
    CO2C2022 = params[:CO2C2022]
    CO2CEACY = params[:CO2CEACY]
    CO2RHFC = params[:CO2RHFC]
    CRDRA = params[:CRDRA]
    CTF = params[:CTF]
    CYRA = params[:CYRA]
    DRC = params[:DRC]
    ECRA22 = params[:ECRA22]
    EGB22 = params[:EGB22]
    @parameters EROCFSP = params[:EROCFSP] [description = "Extra ROC in Food Sector Productivity from 2022 1/y"]
    FCG = params[:FCG]
    FF80 = inits[:CRLA] * params[:FU80]
    FFLREOGRRM = params[:FFLREOGRRM]
    FU80 = params[:FU80]
    FUELER = params[:FUELER]
    FUESQ = params[:FUESQ]
    @parameters GCWR = params[:GCWR] [description = "Goal for Crop Waste Reduction"]
    @parameters GFNRM = params[:GFNRM] [description = "Goal for Fraction New Red Meat"]
    @parameters GFRA = params[:GFRA] [description = "Goal for fraction regenerative agriculture"]
    KCKRM = params[:KCKRM]
    LER80 = params[:LER80]
    MFAM = params[:MFAM]
    OGRR80 = params[:OGRR80]
    OW2022 = params[:OW2022]
    OWEACY = params[:OWEACY]
    ROCFP = params[:ROCFP]
    ROCFSP = params[:ROCFSP]
    SFU = params[:SFU]
    SSP2LMA = params[:SSP2LMA]
    TCTB = params[:TCTB]
    TFFLR = params[:TFFLR]
    UDT = params[:UDT]
    ULP = params[:ULP]

    @variables ACY(t) [description = "Average Crop Yield t-crop/ha/y"]
    @variables ALFL(t) [description = "Acceptable Loss of Forestry Land (1)"]
    @variables AFSRA(t) [description = "Amount of Fertilizer Saved in Reg Ag kgN/ha/y"]
    @variables BALA(t) = inits[:BALA] [description = "BArren LAnd Mha"]
    @variables BIUS(t) [description = "BIofuels USe Mtoe/y"]
    @variables CEM(t) [description = "Cropland Expansion Multiplier (1)"]
    @variables CIRA(t) [description = "Cost Index for Regenerative Agriculture (1)"]
    @variables CO2AFL(t) [description = "CO2 Absorption in Forestry Land GtCO2/y"]
    @variables CO2AFLH(t) [description = "CO2 Absorption in Forestry Land tCO2/ha/y"]
    @variables CO2ELULUC(t) [description = "CO2 Emissions from LULUC GtCO2/y"]
    @variables CO2ELY(t) [description = "CO2 Effect on Land Yield (1)"]
    @variables CO2RFC(t) [description = "CO2 Release from Forest Cut GtCO2/y"]
    @variables COFE(t) [description = "COst of FErtilizer Gdollar/y"]
    @variables COFO(t) [description = "COst of FOod Gdollar/y"]
    @variables CRA(t) [description = "Cost of Regenerative Agriculture Gdollar/y"]
    @variables CRBA(t) [description = "CRop BAlance (1)"]
    @variables CRBI(t) [description = "CRops for BIofuel Mt-crop/y"]
    @variables CRDE(t) [description = "CRop DEmand Mt-crop/y"]
    @variables CREX(t) [description = "CRopland EXpansion Mha/y"]
    @variables CREXR(t) [description = "CRopland EXpansion Rate 1/y"]
    @variables CRLA(t) = inits[:CRLA] [description = "CRopLAnd Mha"]
    @variables CRLO(t) [description = "CRopland LOss Mha/y"]
    @variables CRSU(t) [description = "CRop SUpply (after 20 % waste) Mt-crop/y"]
    @variables CRUS(t) [description = "CRop USe Mt/y"]
    @variables CRUSP(t) [description = "CRop USe per Person t-crop/p/y"]
    @variables CSQCA(t) [description = "Change in Soil Quality in Conv Ag t-crop/ha/y/y"]
    @variables CSRA(t) [description = "Crop Supply Reg Ag Mt-crop/y"]
    @variables CWR(t) [description = "Crop Waste Reduction (1)"]
    @variables DCS(t) [description = "Desired Crop Supply Mt-crop/y"]
    @variables DCSCA(t) [description = "Desired Crop Supply Conv Ag Mt-crop/y"]
    @variables DCYCA(t) [description = "Desired Crop Yield in Conv Ag t-crop/ha/y"]
    @variables DRM(t) [description = "Demand for Red Meat Mt-red-meat/y"]
    @variables DRMP(t) [description = "Demand for Red Meat per Person kg-red-meat/p/y"]
    @variables ECFT(t) [description = "Extra Cost of Food Turnaround Gdollar/y"]
    @variables ECFTSGDP(t) [description = "Extra Cost of Food Turnaround as Share of GDP (1)"]
    @variables ECO2ARA(t) [description = "Extra CO2 Absorption in Reg Ag GtCO2/y"]
    @variables ECRA(t) [description = "Extra Cost of Reg Ag dollar/ha/y"]
    @variables FAM(t) [description = "Forest Absorption Multipler (1)"]
    @variables FCR(t) [description = "Fertilizer Cost Reduction Gdollar/y"]
    @variables FEER(t) [description = "Fertilizer Effect on Erosion Rate (1)"]
    @variables FERM(t) [description = "FEed for Red Meat Mt-crop/y"]
    @variables FEUS(t) [description = "FErtilizer USe Mt/y"]
    @variables FFI(t) [description = "Food Footprint Index (1980=1)"]
    @variables FFLR(t) [description = "Fraction Forestry Land Remaining (1)"]
    @variables FFLREOGRR(t) [description = "FFLReoOGRR"]
    @variables FNRM(t) [description = "Fraction New Red Meat (1)"]
    @variables FOFO(t) [description = "FOod FOotprint"]
    @variables FOLA(t) = params[:FOLA1980] [description = "FOrestry LAnd Mha"]
    @variables FPI(t) [description = "Fertilizer Productivity Index (1980=1)"]
    @variables FRA(t) [description = "Fraction Regenerative Agriculture (1)"]
    @variables FSPI(t) [description = "Food Sector Productivity Index (1980=1)"]
    @variables FUCA(t) [description = "Fertilizer Use in Conv Ag kgN/ha/y"]
    @variables FUP(t) [description = "FErtilizer USe per Person kg/p/y"]
    @variables GLY(t) [description = "Grazing Land Yield kg-red-meat/ha/y"]
    @variables GLY80(t) [description = "Grazing Land Yied in 1980 kg-red-meat/ha/y"]
    @variables GRLA(t) = inits[:GRLA] [description = "GRazing LAnd Mha"]
    @variables IUL(t) [description = "Indicated Urban Land Mha"]
    @variables LERM(t) [description = "Land ERosion Multiplier (1)"]
    @variables LER(t) [description = "Land Erosion Rate 1/y"]
    @variables LFL(t) [description = "Loss of Forest Land Mha/y"]
    @variables LOCR(t) [description = "LOss of CRopland Mha/y"]
    @variables NDRA(t) [description = "Number of Doublings in Reg Ag (1)"]
    @variables NFL(t) [description = "New Forestry Land Mha/y"]
    @variables NGL(t) [description = "New Grazing Land Mha/y"]
    @variables OGFA(t) = inits[:OGFA] [description = "Old Growth Forest Area Mha 1"]
    @variables OGRE(t) [description = "Old Growth Removal Mha/y"]
    @variables OGRR(t) [description = "Old Growth Removal Rate 1/y"]
    @variables OGRRM(t) [description = "Old Growth Removal Rate Multiplier (1)"]
    @variables PCB(t) [description = "Perceived Crop Balance (1)"]
    @variables PRMGL(t) [description = "Potential Red Meat from Grazing Land Mt-red-meat/y"]
    @variables RAA(t) [description = "Regenerative Agriculture Area Mha"]
    @variables RMF(t) [description = "Red Meat from Feedlots Mt-red-meat/y"]
    @variables RMGL(t) [description = "Red Meat from Grazing Land Mt-red-meat/y"]
    @variables RMSP(t) [description = "Red meat Supply per Person kg-red-meat/p/y"]
    @variables ROCSQCA(t) [description = "ROC in Soil Quality in Conv Ag 1/y"]
    @variables SQICA(t) = inits[:SQICA] [description = "Soil Quality Index in Conv Ag (1980=1)"]
    @variables TFA(t) [description = "Total Forest Area Mha"]
    @variables TFUCA(t) [description = "Traditional Fertilizer Use in Conv Ag kgN/ha/y"]
    @variables TUC(t) [description = "Traditional Use of Crops Mt/y"]
    @variables TUCERM(t) [description = "Traditional Use of Crops Ex Red Meat Mt/y"]
    @variables TUCERMP(t) [description = "Traditional Use of Crops Ex Red Meat per Person kg-crop/p/y"]
    @variables TUCP(t) [description = "Traditional Use of Crops per Person kg-crop/p/y"]
    @variables TUFRM(t) [description = "Traditional Use of Feed for Red Meat Mt-crop/y"]
    @variables TURMP(t) [description = "Traditional Use of Red Meat per Person kg-red-meat/p/y"]
    @variables UREX(t) [description = "URban EXpansion Mha/y"]
    @variables URLA(t) = inits[:URLA] [description = "URban LAnd Mha"]
    @variables WELY(t) [description = "Warming Effect on Land Yield (1)"]

    # inventory
    DAT = params[:DAT]
    DDI1980 = params[:DDI1980]
    DIC = params[:DIC]
    DRI = params[:DRI]
    ICPT = params[:ICPT]
    MRIWI = params[:MRIWI]
    OO = params[:OO]
    PH = params[:PH]
    PPU = params[:PPU]
    SAT = params[:SAT]
    INVEODDI = params[:INVEODDI]
    INVEOIN = params[:INVEOIN]
    INVEOSWI = params[:INVEOSWI]
    SRI = params[:SRI]
    SWI1980 = params[:SWI1980]
    TAS = params[:TAS]

    @variables DELDI(t) = inits[:DELDI] [description = "DELivery Delay - Index (1)"]
    @variables EPP(t) = inits[:EPP] [description = "Effective Purchasing Power Gdollar/y"]
    @variables INV(t) = inits[:INV] [description = "INVentory Gu"]
    @variables PRIN(t) = inits[:PRIN] [description = "PRice INdex (1980=1)"]
    @variables PRI(t) = inits[:PRI] [description = "Perceived Relative Inventory (1)"]
    @variables RS(t) = inits[:RS] [description = "Recent Sales Gu/y"]
    @variables SSWI(t) = inits[:SSWI] [description = "ShiftS Worked - Index (1)"]
    @variables CDDI(t) [description = "Change in DDI 1/y"]
    @variables CPI(t) [description = "Change in Price Index 1/y"]
    @variables DEL(t) [description = "DELiveries Gu/y"]
    @variables DSWI(t) [description = "Desired Shifts Worked - Index (1)"]
    @variables GDP(t) [description = "GDP Gdollar/y"]
    @variables IC(t) [description = "Inventory Coverage y"]
    @variables IR(t) [description = "Inflation Rate 1/y"]
    @variables NI(t) [description = "National Income Gdollar/y"]
    @variables OUTP(t) [description = "OUTPut Gu/y"]
    @variables PNIS(t) [description = "Pink Noise In Sales (1)"]
    @variables ROCDDI(t) [description = "ROC in DDI 1/y"]
    @variables SA(t) [description = "SAles Gdollar/y"]

    # labourmarket
    AUR = params[:AUR]
    FIC = params[:FIC]
    GDPP1980 = params[:GDPP1980]
    GDPPEROCCLRM = params[:GDPPEROCCLRM]
    GENLPR = params[:GENLPR]
    NLPR80 = params[:NLPR80]
    PFTJ = params[:PFTJ]
    PFTJ80 = params[:PFTJ80]
    PRUN = params[:PRUN]
    PUELPR = params[:PUELPR]
    ROCECLR80 = params[:ROCECLR80]
    RWER = params[:RWER]
    TAHW = params[:TAHW]
    TELLM = params[:TELLM]
    TENHW = params[:TENHW]
    TYLD = params[:TYLD]
    WSOECLR = params[:WSOECLR]
    WSOELPR = params[:WSOELPR]

    @variables AGIW(t) [description = "Average Gross Income per Worker kdollar/p/y"]
    @variables AHW(t) [description = "Average Hours Worked kh/y"]
    @variables AHW1980(t) [description = "Average Hours Worked in 1980 kh/y"]
    @variables AVWO(t) [description = "AVailable WOrkforce Mp"]
    @variables CECLR(t) [description = "Change in Embedded CLR kcu/ftj/y"]
    @variables CHWO(t) [description = "CHange in WOrkforce Mp/y"]
    @variables CWRA(t) [description = "Change in Wage RAte dollar/ph/y"]
    @variables CWSO(t) [description = "Change in WSO 1/y"]
    @variables ECLR(t) = inits[:ECLR] [description = "Embedded CLR kcu/ftj"]
    @variables ENLPR2022(t) [description = "Extra Normal LPR from 2022 (1)"]
    @variables GDPPEROCCLR(t) [description = "GDPppeoROCCLR"]
    @variables HFD(t) [description = "Hiring/Firing Delay y"]
    @variables HWMGDPP(t) [description = "Hours Worked Mult from GDPpP (1)"]
    @variables IWEOCLR(t) [description = "Indicated Wage Effect on Optimal CLR (1)"]
    @variables ILPR(t) [description = "Indicated Labour Participation Rate (1)"]
    @variables LAPR(t) [description = "LAbour PRoductivity dollar/ph"]
    @variables LAUS(t) = inits[:LAUS] [description = "LAbour USe Gph/y"]
    @variables LAUS80(t) [description = "LAbour USe in 1980 Gph/y"]
    @variables LPR(t) = inits[:ILPR] [description = "Labour Participation Rate (1)"]
    @variables LTEWSO(t) [description = "Long-Term Erosion of WSO 1/y"]
    @variables NHW(t) = inits[:NHW] [description = "Normal Hours Worked kh/ftj/y"]
    @variables NLPR(t) [description = "Normal LPR (1)"]
    @variables OCLR(t) [description = "Optimal Capital Labour Ratio kcu/ftj"]
    @variables OPWO(t) [description = "OPtimal WOrkforce Mp"]
    @variables PART(t) [description = "Participation (1)"]
    @variables PSW(t) [description = "Perceived Surplus Workforce (1)"]
    @variables PURA(t) = inits[:PURA] [description = "Perceived Unemployment RAte (1)"]
    @variables ROCECLR(t) [description = "ROC in ECLR 1/y"]
    @variables ROCWSO(t) [description = "ROC in WSO - Table 1/y"]
    @variables TCT(t) [description = "Time to Change Tooling y"]
    @variables UNEM(t) [description = "UNEMployed Mp"]
    @variables UR(t) [description = "UNemployment RAte (1)"]
    @variables UPT(t) [description = "Unemployment perception time y"]
    @variables WAP(t) [description = "Working Age Population Mp"]
    @variables WARA(t) = inits[:WARA] [description = "WAge RAte dollar/ph"]
    @variables WASH(t) [description = "WAge Share (1)"]
    @variables WEOCLR(t) = inits[:WEOCLR] [description = "Wage Effect on Optimal CLR (1)"]
    @variables WF(t) = inits[:WF] [description = "WorkForce Mp"]
    @variables WRER(t) [description = "Wage Rate Erosion Rate 1/y"]
    @variables WRE(t) [description = "Wage Rate Erosion dollar/ph/y"]
    @variables WSO(t) = inits[:WSO] [description = "Worker Share of Output (1)"]

    # other
    INELOK = params[:INELOK]
    NK = params[:NK]
    TEGR = params[:TEGR]

    @variables CFETA(t) [description = "Cost of Food and Energy TAs GDollar/y"]
    @variables CTA(t) [description = "Cost of TAs GDollar/y"]
    @variables FB15(t) [description = "Fraction Below 15 kDollar/p/y (1)"]
    @variables IEL(t) [description = "Inequity Effect on Logistic k (1)"]
    @variables LK(t) [description = "Logistic K (1)"]
    @variables PB15(t) [description = "Population Below 15 kDollar/p/y Mp"]
    @variables PGDPP(t) = inits[:PGDPP] [description = "Past GDP per Person kDollar/y"]
    @variables RGGDPP(t) [description = "Rate of Growth in GDP per Person 1/y"]

    # output
    CAPPIS1980 = params[:CAPPIS1980]
    CAPPUS1980 = params[:CAPPUS1980]
    CC1980 = params[:CC1980]
    CTPIS = params[:CTPIS]
    CTPUS = params[:CTPUS]
    EMCUC = params[:EMCUC]
    FCI = params[:FCI]
    LCPIS1980 = params[:LCPIS1980]
    OW2022 = params[:OW2022]
    OWECCM = params[:OWECCM]
    OWELCM = params[:OWELCM]
    @parameters USPIS2022 = params[:USPIS2022] [description = "Unconventional Stimulus in PIS from 2022 (share of GDP)"]
    @parameters USPUS2022 = params[:USPUS2022] [description = "Unconventional Stimulus in PUS from 2022 (share of GDP)"]
    CBCEFRA = params[:CBCEFRA]
    CUCPIS1980 = (params[:CAPPIS1980] / params[:LCPIS1980]) * params[:CTPIS] * params[:EMCUC]
    ED1980 = params[:ED1980]
    EDEFRA = params[:EDEFRA]
    EDELCM = params[:EDELCM]
    FRA1980 = params[:FRA1980]
    FRACAM = params[:FRACAM]
    GDPP1980 = params[:GDPP1980]
    GDPPEFRACA = params[:GDPPEFRACA]
    IPT = params[:IPT]
    JOBS1980 = params[:JOBS1980]
    KAPPA = params[:KAPPA]
    LAMBDA = params[:LAMBDA]
    LAUS1980 = params[:LAUS1980]
    OG1980 = params[:OG1980]
    OO1980 = params[:OO1980]
    PRUN = params[:CU1980] * (1 + params[:MA1980])
    TOED = params[:TOED]
    WSOEFRA = params[:WSOEFRA]

    @variables AVCA(t) [description = "AVailable CApital Gdollar/y"]
    @variables CAPIS(t) [description = "Capacity Addition PIS Gcu/y"]
    @variables CAPUS(t) [description = "Capacity Addition PUS Gcu/y"]
    @variables CDPIS(t) [description = "Capacity Discard PIS Gcu/y"]
    @variables CDPUS(t) [description = "Capacity Discard PUS Gcu/y"]
    @variables CIPIS(t) [description = "Capacity Initiation PIS Gcu/y"]
    @variables CIPUS(t) [description = "Capacity Initiation PUS Gcu/y"]
    @variables COCA(t) [description = "COst of CApacity dollar/cu"]
    @variables CPIS(t) = params[:CAPPIS1980] [description = "Capacity PIS Gcu"]
    @variables CPUS(t) = params[:CAPPUS1980] [description = "Capacity PUS Gcu"]
    @variables CRR(t) [description = "Capacity Renewal Rate 1/y"]
    @variables CUCPIS(t) = inits[:CUCPIS] [description = "Capacity Under Construction PIS Gcu"]
    @variables CUCPIS1980(t) [description = "CUC PIS in 1980 Gcu"]
    @variables CUCPUS(t) = inits[:CUCPUS] [description = "Capacity Under Construction PUS Gcu"]
    @variables CUCPUS1980(t) [description = "CUC PUS in 1980 Gcu"]
    @variables ECR(t) [description = "Effect of Capacity Renewal 1/y"]
    @variables ETFP(t) = inits[:ETFP] [description = "Embedded TFP (1)"]
    @variables INCPIS(t) [description = "Investment in New Capacity PIS Gdollar/y"]
    @variables LCPUS(t) [description = "Life of Capacity PUS y"]
    @variables LCPUS1980(t) [description = "Life of Capacity PUS in 1980 y"]
    @variables OBSGIPIS(t) [description = "Off-Balance Sheet Govmnt Inv in PIS (share of GDP)"]
    @variables OBSGIPUS(t) [description = "Off-Balance-Sheet Govmnt Inv in PUS (share of GDP)"]
    @variables OWECC(t) [description = "OWeoCOC (1)"]
    @variables OWELC(t) [description = "OWeoLOC (1)"]
    @variables CBCEFCA(t) [description = "CBC Effect on Flow to Capacity Addion (1)"]
    @variables CAP(t) [description = "CAPacity Gcu"]
    @variables EDE(t) [description = "Excess DEmand (1)"]
    @variables EDEFCA(t) [description = "ED Effect on Flow to Capacity Addition (1)"]
    @variables EDELC(t) [description = "Excess Demand Effect on Life of Capacity (1)"]
    @variables FACNC(t) = inits[:FACNC] [description = "Fraction of Available Capital to New Capacity (1)"]
    @variables FRACAMGDPPL(t) [description = "FRACA Mult from GDPpP - Line (1)"]
    @variables ISGDP(t) [description = "Investment Share of GDP (1)"]
    @variables LCPIS(t) [description = "Life of Capacity PIS y"]
    @variables OGR(t) = params[:OG1980] [description = "Output Growth Rate 1/y"]
    @variables OLY(t) = inits[:OLY] [description = "Output Last Year Gdollar/y"]
    @variables OOV(t) [description = "Optimal Ouput - Value Gdollar/y"]
    @variables ORO(t) [description = "Optimal Real Output Gu/y"]
    @variables PEDE(t) = 1 [description = "Perceived Excess DEmand (1)"]
    @variables WSOEFCA(t) [description = "WSO Effect on Flow to Capacity Addition (1)"]

    # population
    CMFR = params[:CMFR]
    DNC80 = params[:DNC80]
    DNCA = params[:DNCA]
    DNCG = params[:DNCG]
    DNCM = params[:DNCM]
    EIP = params[:EIP]
    EIPF = params[:EIPF]
    EPA22 = params[:EPA22]
    FADFS = params[:FADFS]
    FP = params[:FP]
    FW = params[:FW]
    GDPP1980 = params[:GDPP1980]
    @parameters GEFR = params[:GEFR] [description = "Goal for Extra Fertility Reduction"]
    GEPA = params[:GEPA]
    LE1980 = params[:LE1980]
    LEA = params[:LEA]
    LEEPA = params[:LEEPA]
    LEG = params[:LEG]
    LEMAX = params[:LEMAX]
    MFM = params[:MFM]
    MLEM = params[:MLEM]
    OW2022 = params[:OW2022]
    OWELE = params[:OWELE]
    PA1980 = params[:PA1980]
    SSP2FA2022F = params[:SSP2FA2022F]
    TAHI = params[:TAHI]

    @variables A0020(t) = inits[:A0020] [description = "Aged 0-20 years Mp"]
    @variables A2040(t) = inits[:A2040] [description = "Aged 20-40 years Mp"]
    @variables A20PA(t) [description = "Aged 20-Pension Age Mp"]
    @variables A4060(t) = inits[:A4060] [description = "Aged 40-60 Mp"]
    @variables A60PL(t) = inits[:A60PL] [description = "Aged 60 + Mp"]
    @variables BIRTHR(t) [description = "Birth Rate 1/y"]
    @variables BIRTHS(t) [description = "Births Mp/y"]
    @variables CEFR(t) [description = "Cost of Extra Fertility Reduction (share of GDP)"]
    @variables DEATHR(t) [description = "Death Rate 1/y"]
    @variables DEATHS(t) [description = "Deaths Mp/y"]
    @variables DNC(t) [description = "Desired No of Children 1"]
    @variables DR(t) [description = "Dependency Ratio p/p"]
    @variables EFR(t) [description = "Extra Fertility Reduction (1)"]
    @variables EGDPP(t) = inits[:EGDPP] [description = "Effective GDP per Person kDollar/p/y"]
    @variables EPA(t) = inits[:EPA] [description = "Extra Pension Age y"]
    @variables FM(t) [description = "Fertility Multiplier (1)"]
    @variables GDPP(t) [description = "GDP per Person kDollar/p/y"]
    @variables LE(t) = inits[:LE] [description = "Life Expectancy y"]
    @variables LE60(t) [description = "LE at 60 y"]
    @variables LEM(t) [description = "Life Expectancy Multipler (1)"]
    @variables (LV_DEATHS(t))[1:ORDER] = fill(inits[:DEATHS] * (inits[:LE] - 60) / ORDER, ORDER) [description = "LV functions for deaths Mp/y"]
    @variables (LV_PASS20(t))[1:ORDER] = fill(inits[:PASS20] * 20 / ORDER, ORDER) [description = "LV functions for passing 20 Mp/y"]
    @variables (LV_PASS40(t))[1:ORDER] = fill(inits[:PASS40] * 20 / ORDER, ORDER) [description = "LV functions for passing 40 Mp/y"]
    @variables (LV_PASS60(t))[1:ORDER] = fill(inits[:PASS60] * 20 / ORDER, ORDER) [description = "LV functions for passing 60 Mp/y"]
    @variables OF(t) [description = "Observed Fertility 1"]
    @variables OP(t) [description = "On Pension Mp"]
    @variables PA(t) = inits[:PA] [description = "Pension Age y"]
    @variables PASS20(t) [description = "Passing 20 Mp/y"]
    @variables PASS40(t) = inits[:PASS40] [description = "Passing 40 Mp/y"]
    @variables PASS60(t) = inits[:PASS60] [description = "Passing 60 Mp/y"]
    @variables PGR(t) [description = "Population Growth Rate 1/y"]
    @variables POP(t) [description = "Population Mp"]
    @variables PW(t) [description = "Pensioners per Worker p/p"]
    @variables (RT_DEATHS(t))[1:ORDER] = fill(inits[:DEATHS] * (inits[:LE] - 60) / ORDER, ORDER) [description = "RT functions for deaths Mp/y"]
    @variables (RT_PASS20(t))[1:ORDER] = fill(inits[:PASS20] * 20 / ORDER, ORDER) [description = "RT functions for passing 20 Mp/y"]
    @variables (RT_PASS40(t))[1:ORDER] = fill(inits[:PASS40] * 20 / ORDER, ORDER) [description = "RT functions for passing 40 Mp/y"]
    @variables (RT_PASS60(t))[1:ORDER] = fill(inits[:PASS60] * 20 / ORDER, ORDER) [description = "RT functions for passing 60 Mp/y"]
    @variables WELE(t) [description = "Warming Effect on Life Expectancy (1)"]

    # public
    OW2022 = params[:OW2022]
    CTA2022 = params[:CTA2022]
    CTPIS = params[:CTPIS]
    DROTA1980 = params[:DROTA1980]
    EDROTA2022 = params[:EDROTA2022]
    FUATA = params[:FUATA]
    GDPTL = params[:GDPTL]
    IIEEROTA = params[:IIEEROTA]
    IPR1980 = params[:IPR1980]
    IPRVPSS = params[:IPRVPSS]
    IPT = params[:IPT]
    @parameters MIROTA2022 = params[:MIROTA2022] [description = "Max Imported ROTA from 2022 1/y"]
    OWETFP = params[:OWETFP]
    SC1980 = params[:SC1980]
    SCROTA = params[:SCROTA]
    XETAC2022 = params[:XETAC2022]
    XETAC2100 = params[:XETAC2100]

    @variables CTFP(t) [description = "Change in TFP 1/y"]
    @variables DRTA(t) [description = "Domestic Rate of Technological Advance 1/y"]
    @variables ECTAF2022(t) [description = "Extra Cost of TAs From 2022 GDollar/y"]
    @variables ECTAGDP(t) [description = "Extra Cost of TAs as share of GDP (1)"]
    @variables GSSGDP(t) [description = "Govmnt Spending as Share of GDP"]
    @variables IPR(t) [description = "Infrastructure Purchases Ratio y"]
    @variables IROTA(t) [description = "Imported ROTA 1/y"]
    @variables ITFP(t) [description = "Indicated TFP (1)"]
    @variables OWTFP(t) [description = "OWeoTFP"]
    @variables PLUA(t) [description = "Productivity Loss from Unprofitable Activity (1)"]
    @variables PPP(t) [description = "Productivity of Public Purchases (1)"]
    @variables PSEP(t) [description = "Public SErvices per Person kdollar/p/y"]
    @variables PSP(t) [description = "Public Spending per Person kDollar/p/y"]
    @variables RROTAI(t) [description = "Reduction in ROTA from Inequality 1/y"]
    @variables RTA(t) [description = "Rate of Technological Advance 1/y"]
    @variables RTFPUA(t) = inits[:RTFPUA] [description = "Reduction in TFP from Unprofitable Activity (1)"]
    @variables SC(t) [description = "State Capacity (fraction of GDP)"]
    @variables TFPEE5TA(t) = inits[:TFPEE5TA] [description = "TFP Excluding Effect of 5TAs (1)"]
    @variables TFPIE5TA(t) [description = "TFP Including Effect of 5TAs (1)"]
    @variables VPSS(t) [description = "Value of Public Services Supplied GDollar/y"]

    # wellbeing
    AI = params[:AI]
    AP = params[:AP]
    AWBPD = params[:AWBPD]
    DRDI = params[:DRDI]
    DRPS = params[:DRPS]
    EIP = params[:EIP]
    EIPF = params[:EIPF]
    GWEAWBGWF = params[:GWEAWBGWF]
    IEAWBIF = params[:IEAWBIF]
    MWBGW = params[:MWBGW]
    NRD = params[:NRD]
    PAEAWBF = params[:PAEAWBF]
    PREAWBF = params[:PREAWBF]
    PESTF = params[:PESTF]
    SPS = params[:SPS]
    STEERDF = params[:STEERDF]
    STRERDF = params[:STRERDF]
    TCRD = params[:TCRD]
    TDI = params[:TDI]
    TEST = params[:TEST]
    TI = params[:TI]
    THPA = params[:THPA]
    TPR = params[:TPR]
    TPS = params[:TPS]
    TW = params[:TW]

    @variables AWBDI(t) [description = "Average WellBeing from Disposable Income (1)"]
    @variables AWBGW(t) [description = "Average WellBeing from Global Warming (1)"]
    @variables AWBI(t) [description = "Average WellBeing Index (1)"]
    @variables AWBIN(t) [description = "Average WellBeing from INequality (1)"]
    @variables AWBPS(t) [description = "Average wellBeing from public spending (1)"]
    @variables AWBP(t) [description = "Average WellBeing from Progress (1)"]
    @variables IEST(t) [description = "Inequity Effect on Social Trust (1)"]
    @variables IPP(t) [description = "Introduction Period for Policy y"]
    @variables IRD(t) [description = "Indicated Reform Delay y"]
    @variables IST(t) [description = "Indicated Social Trust (1)"]
    @variables ORP(t) = inits[:ORP] [description = "Observed Rate of Progress 1/y"]
    @variables PAWBI(t) = inits[:PAWBI] [description = "Past AWI (1)"]
    @variables PSESTR(t) [description = "Public Spending Effect on Social TRust (1)"]
    @variables PSSGDP(t) [description = "Public Spending as Share of GDP"]
    @variables RD(t) = inits[:RD] [description = "Reform Delay y"]
    @variables SOTR(t) = inits[:SOTR] [description = "Social Trust (1)"]
    @variables STE(t) = inits[:STE] [description = "Social TEnsion (1)"]
    @variables STEERD(t) [description = "Social TEnsion Effect on Reform Delay (1)"]
    @variables STRERD(t) [description = "Social TRust Effect on Reform Delay (1)"]
    @variables WBEP(t) [description = "WellBeing Effect of Participation (1)"]

    #other parameters 
    FOLA1980 = params[:FOLA1980]
    NHW1980 = params[:NHW1980]
    WSO1980 = params[:WSO1980]
    WF1980 = params[:WF1980]
    STE1980 = params[:STE1980]
    STR1980 = params[:STR1980]


    # climate
    add_equation!(eqs, KN2OEKF ~ KN2OKF1980 * exp(-(RDN2OKF) * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-(ERDN2OKF2022) * (t - 2022)), 1))
    add_equation!(eqs, MMN2OE ~ FEUS * KN2OEKF / 1000)
    add_equation!(eqs, NN2OE ~ withlookup(t, [(1980.0, 0.009), (2020.0, 0.009), (2099.27, 0.0)]))
    add_equation!(eqs, N2OE ~ NN2OE + MMN2OE)
    add_equation!(eqs, D(N2OA) ~ N2OE - N2OBD)
    add_equation!(eqs, N2OBD ~ N2OA / LN2OA)
    add_equation!(eqs, N2OCA ~ N2OA / GN2OPP)
    add_equation!(eqs, N2OFPP ~ withlookup(t, [(1980.0, 0.43), (2000.0, 0.64), (2010.0, 0.73), (2020.0, 0.8), (2100.0, 1.0)]))
    add_equation!(eqs, FN2O ~ N2OCA * N2OFPP)
    add_equation!(eqs, KCH4EKC ~ KCH4KC1980 * exp(-(RDCH4KC) * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-(ERDCH4KC2022) * (t - 2022)), 1))
    add_equation!(eqs, MMCH4E ~ CRSU * KCH4EKC / 1000)
    add_equation!(eqs, NCH4E ~ withlookup(t, [(1980.0, 0.19), (2020.0, 0.19), (2100.0, 0.19)]))
    add_equation!(eqs, CH4E ~ NCH4E + MMCH4E)
    add_equation!(eqs, D(CH4A) ~ CH4E - CH4BD)
    add_equation!(eqs, CH4BD ~ CH4A / LCH4A)
    add_equation!(eqs, CH4CA ~ CH4A / GCH4PP)
    add_equation!(eqs, CH4FPP ~ withlookup(t, [(1980.0, 0.82), (2000.0, 0.94), (2020.0, 1.01), (2100.0, 1.1)]))
    add_equation!(eqs, FCH4 ~ CH4CA * CH4FPP)
    add_equation!(eqs, OWLCO2 ~ IfElse.ifelse(t > 2022, 1 + SOWLCO2 * (OW / OBWA2022 - 1), 1))
    add_equation!(eqs, LECO2A ~ LECO2A1980 * OWLCO2)
    add_equation!(eqs, CO2FCH4 ~ CH4BD * TCO2PTCH4)
    add_equation!(eqs, CO2AB ~ (CO2A - CO2A1850) / LECO2A)
    add_equation!(eqs, CO2E ~ CO2EI + CO2ELULUC - DACCO2)
    add_equation!(eqs, CO2GDP ~ (CO2E / GDP) * 1000)
    add_equation!(eqs, CAC ~ DACCO2 * CCCSt)
    add_equation!(eqs, DACCO2 ~ IfElse.ifelse(t > 2022, ramp(t, (DACCO22100) / IPP, 2022, 2022 + IPP), 0))
    add_equation!(eqs, D(CO2A) ~ CO2E - CO2AB + 2 * CO2FCH4) #STRANGE EQUATION!
    add_equation!(eqs, CO2CA ~ CO2A / GCO2PP)
    add_equation!(eqs, CO2FPP ~ withlookup(t, [(1980.0, 0.0032), (1990.0, 0.0041), (2000.0, 0.0046), (2020.0, 0.0051), (2100.0, 0.006)]))
    add_equation!(eqs, FCO2 ~ CO2CA * CO2FPP)
    add_equation!(eqs, FOG ~ withlookup(t, [(1980.0, 0.18), (2000.0, 0.36), (2020.0, 0.39), (2050.0, 0.37), (2100.0, 0.0)]))
    add_equation!(eqs, MMF ~ FCO2 + FOG + FCH4 + FN2O)
    add_equation!(eqs, GHGE ~ CO2E * TCO2ETCO2 + CH4E * TCO2ETCH4 + N2OE * TCO2ETN2O)
    add_equation!(eqs, AL1980 ~ (ISCEGA1980 * ALIS + (GLSU - ISCEGA1980) * ALGAV) / GLSU)
    add_equation!(eqs, AL ~ (ISCEGA * ALIS + (GLSU - ISCEGA) * ALGAV) / GLSU)
    add_equation!(eqs, TRHGS ~ (TRSS1980 * ((OW + 297) / 297)) * (AL / AL1980))
    add_equation!(eqs, HTS ~ EHS * TRHGS) # EHS
    add_equation!(eqs, MRS ~ MRS1980 * (OW / WA1980))
    add_equation!(eqs, MRDI ~ MRS / SVDR)
    add_equation!(eqs, ECIM ~ MRDI * AI1980 * TPM3I * HRMI)
    add_equation!(eqs, MEL ~ ISCEGA * MRS)
    add_equation!(eqs, D(ISCEGA) ~ -MEL)
    add_equation!(eqs, ISC ~ ISCEGA * 100)
    add_equation!(eqs, WVC ~ WVC1980 * (1 + OWWV * (OW / WA1980 - 1)))
    add_equation!(eqs, WVF ~ WVF1980 * (1 + WVWVF * (WVC / WVC1980 - 1)))
    add_equation!(eqs, TMMF ~ MMF + WVF)
    add_equation!(eqs, EWFF ~ (TMMF * GLSU) * 31.5 / 1000)
    add_equation!(eqs, OW ~ WA1980 + (EHS - EH1980) * WFEH)
    add_equation!(eqs, REHE ~ withlookup(OW, [(0.0, 1.0), (1.2, 4.8), (2.0, 8.6), (2.9, 14.0), (5.2, 40.0)]))
    add_equation!(eqs, TRHGA ~ TRSA1980 * ((OW + 287) / 287))
    add_equation!(eqs, HDO ~ EHS * TRHGA)
    add_equation!(eqs, D(EHS) ~ EWFF - ECIM - HDO - HTS)
    smooth!(eqs, PWA, OW, PD)

    # demand
    add_equation!(eqs, BCIL ~ CFWB + CFGB)
    add_equation!(eqs, BCISNI ~ BCIL / NI)
    add_equation!(eqs, BITRO ~ min(1, ITRO1980) + ramp(t, (ITRO2022 - ITRO1980) / 42, 1980, 2022) + ramp(t, (GITRO - ITRO2022) / 78, 2022, 2100))
    add_equation!(eqs, CANCD ~ pulse(t, 2022, 1) * GD * FGDC2022)
    add_equation!(eqs, CD ~ WCD - STW + OC - STO)
    add_equation!(eqs, CFGB ~ GIC + GP - GND)
    add_equation!(eqs, CFWB ~ WIC + WP - WND)
    add_equation!(eqs, CPP ~ CD / POP)
    add_equation!(eqs, CSGDP ~ CD / NI)
    add_equation!(eqs, CONTR ~ CSGDP + GSGDP + SSGDP)
    add_equation!(eqs, EGTF2022 ~ IfElse.ifelse(t > 2022, EGTRF2022 + EETF2022 + EPTF2022, 0) * NI)
    smooth!(eqs, ETF2022, GETF2022, TINT)
    add_equation!(eqs, ETTAF2022 ~ IfElse.ifelse(t > 2022, ECTAF2022 * FETACPET, 0))
    smooth!(eqs, FGBW, GFGBW, TINT)
    add_equation!(eqs, GCIN ~ GNI - CFGB)
    add_equation!(eqs, D(GD) ~ GND - CANCD - GP)
    add_equation!(eqs, GDB ~ GD / NI)
    add_equation!(eqs, GETF2022 ~ EGTF2022 + ETTAF2022)
    add_equation!(eqs, GFGBW ~ FT1980 + IfElse.ifelse(t > 2022, ETGBW, 0))
    add_equation!(eqs, GFSNI ~ (GIC + GP) / NI)
    add_equation!(eqs, GGI ~ WT + OT + STO + STW + IC2022)
    add_equation!(eqs, GGIS ~ GGI / NI)
    add_equation!(eqs, GIC ~ GD * GBC)
    add_equation!(eqs, GIPC ~ PGCIN - GPU)
    add_equation!(eqs, GND ~ max(0, (MGD - GD) / GDDP) + step(t, GSF2022, 2022) * NI)
    add_equation!(eqs, GNI ~ (WT + OT + STO + STW + IC2022) - TP + ST)
    add_equation!(eqs, GNISNI ~ GNI / NI)
    add_equation!(eqs, GP ~ GD / GPP)
    add_equation!(eqs, GPU ~ PGCIN * GCF)
    add_equation!(eqs, GS ~ GPU + GIPC)
    add_equation!(eqs, GSGDP ~ GS / NI)
    add_equation!(eqs, IC2022 ~ NI * IfElse.ifelse(t > 2022, ramp(t, GEIC / IPP, 2022, 2020 + IPP), 0))
    add_equation!(eqs, INEQ ~ OOIAT / WIAT)
    add_equation!(eqs, INEQI ~ INEQ / INEQ1980)
    add_equation!(eqs, ITO ~ BITRO * NI * (1 - WSO))
    add_equation!(eqs, ITW ~ BITRW * NI * WSO)
    add_equation!(eqs, MGD ~ NI * MGDB)
    add_equation!(eqs, MWD ~ WI * MWDB)
    add_equation!(eqs, OC ~ POCI * OCF)
    add_equation!(eqs, OCF ~ 1 - OSF)
    add_equation!(eqs, OCIN ~ OOIAT)
    add_equation!(eqs, OI ~ NI * (1 - WSO))
    add_equation!(eqs, OOIAT ~ OI - OT)
    add_equation!(eqs, OS ~ POCI - OC)
    add_equation!(eqs, OSF ~ OSF1980 * (1 + GDPOSR * (EGDPP / GDPP1980 - 1)))
    add_equation!(eqs, OT ~ ITO + ETF2022 * FETPO)
    add_equation!(eqs, OTR ~ OT / OI)
    smooth!(eqs, PGCIN, GCIN, TAB)
    smooth!(eqs, POCI, OCIN, TAOC)
    smooth!(eqs, PWCIN, WCIN, TAWC)
    add_equation!(eqs, SSGDP ~ TS / NI)
    add_equation!(eqs, ST ~ STW + STO)
    add_equation!(eqs, STO ~ OC * STR)
    add_equation!(eqs, STW ~ WCD * STR)
    add_equation!(eqs, TP ~ (WT + OT + STO + STW + IC2022) * FGBW)
    add_equation!(eqs, TPP ~ WCIN + GCIN + OCIN - ST)
    add_equation!(eqs, TS ~ OS + WS)
    add_equation!(eqs, WCD ~ PWCIN * WCF)
    add_equation!(eqs, WCIN ~ WIAT - CFWB)
    add_equation!(eqs, D(WD) ~ WND - WP)
    add_equation!(eqs, WDI ~ PWCIN / WF)
    add_equation!(eqs, WFCSI ~ CFWB / WIAT)
    add_equation!(eqs, WI ~ NI * WSO)
    add_equation!(eqs, WIAT ~ WI - WT + TP)
    add_equation!(eqs, WIC ~ WD * WBC)
    add_equation!(eqs, WT ~ ITW + ETF2022 * (1 - FETPO))
    add_equation!(eqs, WTR ~ WT / WI)
    add_equation!(eqs, WND ~ max(0, (MWD - WD) / WDP))
    add_equation!(eqs, WP ~ WD / WPP)
    add_equation!(eqs, WDB ~ WD / WIAT)
    add_equation!(eqs, WS ~ PWCIN - WCD)

    # energy
    add_equation!(eqs, CNEL ~ NEP * CNED)
    add_equation!(eqs, NFCO2PP ~ MNFCO2PP * (1 - exp(-(GDPP / 10))))
    add_equation!(eqs, CO2NFIP ~ (NFCO2PP / 1000) * POP * (1 - FCO2SCCS))
    add_equation!(eqs, FCO2SCCS ~ FCO2SCCS2022 + ramp(t, (GFCO2SCCS - FCO2SCCS2022) / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, CO2EI ~ CO2EP + CO2NFIP)
    add_equation!(eqs, CO2EP ~ UFF * (TCO2PT / 1000) * (1 - FCO2SCCS))
    add_equation!(eqs, CO2EMPP ~ (CO2EI / POP) * 1000)
    add_equation!(eqs, CCCSG ~ CCCSt * ICCSC)
    add_equation!(eqs, ICCSC ~ FCO2SCCS * (CO2NFIP + CO2EP) / (1 - FCO2SCCS))
    add_equation!(eqs, TCO2PT ~ 2.8 * exp(ROCTCO2PT * (t - 1980)))
    add_equation!(eqs, D(EEPI2022) ~ IEEPI)
    add_equation!(eqs, IEEPI ~ EROCEPA2022 * 0 + step(t, EROCEPA2022, 2022))
    add_equation!(eqs, TPPUEBEE ~ withlookup(GDPP, [(0.0, 0.0), (10.0, 4.0), (20.0, 7.0), (30.0, 9.0), (50.0, 12.0), (65.0, 13.0)]))
    add_equation!(eqs, TPPUFFNEUBEE ~ withlookup(GDPP, [(0.0, 0.3), (15.0, 2.0), (25.0, 3.1), (35.0, 4.0), (50.0, 5.0)]))
    add_equation!(eqs, DEBNE ~ (POP * TPPUEBEE * exp(-NIEE * (t - 1980))) / EEPI2022)
    add_equation!(eqs, DFFNEUBNE ~ (POP * TPPUFFNEUBEE * exp(-NIEE * (t - 1980))) / EEPI2022)
    add_equation!(eqs, FNE ~ FNE1980 + ramp(t, (FNE2022 - FNE1980) / 42, 1980, 2022) + ramp(t, (GFNE - FNE2022) / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, ERDNEFFFNE ~ FNE * DFFNEUBNE)
    add_equation!(eqs, CNE ~ (ECRUNEFF / 1000) * ERDNEFFFNE)
    add_equation!(eqs, EIDEFNE ~ ERDNEFFFNE * EUEPRUNEFF)
    add_equation!(eqs, DFFFNEU ~ DFFNEUBNE - ERDNEFFFNE)
    add_equation!(eqs, UFF ~ DFFFNEU + FFE)
    add_equation!(eqs, DE ~ DEBNE + EIDEFNE)
    add_equation!(eqs, DRES ~ REFF1980 + ramp(t, (REFF2022 - REFF1980) / 42, 1980, 2022) + ramp(t, (GREF - REFF2022) / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, DSRE ~ DE * DRES)
    add_equation!(eqs, DREC ~ DSRE / RCUT)
    add_equation!(eqs, DRECC ~ DREC - REC)
    add_equation!(eqs, D(REC) ~ AREC - DIREC)
    add_equation!(eqs, AREC ~ max(0, (DRECC / RECT) + (DIREC)))
    add_equation!(eqs, DIREC ~ REC / LREC)
    add_equation!(eqs, ASWC ~ AREC)
    add_equation!(eqs, D(ACSWCF1980) ~ ASWC)
    add_equation!(eqs, NDSWC ~ log(2) + log(ACSWCF1980 / SWC1980))
    add_equation!(eqs, CISWC ~ (1 - CRDSWC)^NDSWC)
    add_equation!(eqs, CAPEXRED ~ CAPEXRE1980 * CISWC)
    add_equation!(eqs, CAPEXREG ~ CAPEXRED * AREC)
    add_equation!(eqs, OPEXREG ~ OPEXRED * REP)
    add_equation!(eqs, CRE ~ CAPEXREG + OPEXREG)
    add_equation!(eqs, CAPEXFEG ~ CAPEXFED * AFEC)
    add_equation!(eqs, OPEXFEG ~ OPEXFED * FEP)
    add_equation!(eqs, CFE ~ CAPEXFEG + OPEXFEG)
    add_equation!(eqs, CEL ~ CFE + CRE + CNEL)
    add_equation!(eqs, REP ~ REC * RCUT)
    add_equation!(eqs, GHMH2 ~ REP * FREH / KWEPKGH2)
    add_equation!(eqs, GHMt ~ GHMH2 * TPTH2)
    add_equation!(eqs, RHP ~ BEM + GHMt)
    add_equation!(eqs, TWEPEJEE ~ TWHPEJCE * EFPP)
    add_equation!(eqs, IIASAREP ~ REP / TWEPEJEE + RHP / MTPEJCE)
    add_equation!(eqs, FTWEPMt ~ TWEPEJEE / MTPEJCE)
    add_equation!(eqs, IIASAFEP ~ UFF / MTPEJCE)
    add_equation!(eqs, LCEP ~ REP + NEP)
    add_equation!(eqs, DFE ~ max(0, DE - LCEP))
    add_equation!(eqs, DFEC ~ DFE / EKHPY)
    add_equation!(eqs, DFECC ~ (DFEC - FEC) / FECCT + DIFEC)
    add_equation!(eqs, AFEC ~ max(0, DFECC))
    add_equation!(eqs, D(FEC) ~ AFEC - DIFEC)
    add_equation!(eqs, LFEC ~ NLFEC * FCUTLOFC)
    add_equation!(eqs, DIFEC ~ FEC / LFEC)
    add_equation!(eqs, FCUT ~ DFE / FEC)
    add_equation!(eqs, FCUTLOFC ~ 1 + sFCUTLOFC * ((FCUT / EKHPY) - 1))
    add_equation!(eqs, FEP ~ FEC * FCUT)
    add_equation!(eqs, NC ~ withlookup(t, [(1980.0, 75.0), (2000.0, 310.0), (2020.0, 310.0), (2098.9, 310.0)]))
    add_equation!(eqs, NEP ~ NC * NCUT)
    add_equation!(eqs, EP ~ FEP + NEP + REP)
    add_equation!(eqs, ELB ~ EP / DE)
    add_equation!(eqs, FFPNE ~ (FEP + NEP) / EP)
    add_equation!(eqs, EU ~ DFFFNEU + EP / FTWEPMt + RHP)
    add_equation!(eqs, EUPP ~ EU / POP)
    add_equation!(eqs, FFE ~ FEP / FTWEPMt)
    add_equation!(eqs, TCEG ~ (DEBNE * TCE / 1000) * AFMCM)
    add_equation!(eqs, TCFFFNEUG ~ (DFFNEUBNE * TCFFFNEU / 1000) * AFMCM)
    add_equation!(eqs, CFFFNEU ~ (DFFFNEU * TCFFFNEU) / 1000)
    add_equation!(eqs, CG ~ EP * TC)
    add_equation!(eqs, TGC ~ DEBNE * TC)
    add_equation!(eqs, TCEN ~ TCEG + TCFFFNEUG + TGC)
    add_equation!(eqs, TCENSGDP ~ TCEN / GDP)
    add_equation!(eqs, CE ~ CFFFNEU + CEL + CG + CNE + CCCSG + CAC)
    add_equation!(eqs, RECTEC ~ CE / TCEN)
    add_equation!(eqs, CESGDP ~ CE / GDP)
    add_equation!(eqs, ECETSGDP ~ IfElse.ifelse(t > 2022, (CE - TCEN) / GDP, 0))

    # finance
    add_equation!(eqs, CBC ~ CCSD + NCCR)
    add_equation!(eqs, CBC1980 ~ NSR + NBBM + NBOM + NCCR)
    smooth!(eqs, CCSD, TIR + NBOM, FSRT)
    add_equation!(eqs, D(CBSR) ~ CSR)
    add_equation!(eqs, CSR ~ (ISR - CBSR) / SRAT)
    smooth!(eqs, ELTI, PI, IEFT)
    add_equation!(eqs, GBC ~ TIR)
    add_equation!(eqs, ISR ~ NSR * (1 + INSR * (PI / IT - 1) + UNSR * (PU / UT - 1)))
    add_equation!(eqs, NCCR ~ 0.02 * (1 + GRCR * (OGR / 0.03 - 1)))
    smooth!(eqs, PI, IR, IPTCB)
    smooth!(eqs, PU, UR, UPTCB)
    add_equation!(eqs, TGIR ~ GBC + ELTI)
    add_equation!(eqs, TIR ~ CBSR + NBBM)
    add_equation!(eqs, WBC ~ CCSD)

    # foodland
    add_equation!(eqs, ACY ~ (DCYCA * (1 - FRA) * SQICA + CYRA * FRA) * CO2ELY * WELY)
    add_equation!(eqs, ALFL ~ 1 - exp(-FFLR / TFFLR))
    add_equation!(eqs, AFSRA ~ 268 - SFU)
    add_equation!(eqs, D(BALA) ~ CRLO)
    add_equation!(eqs, BIUS ~ withlookup(t, [(1980.0, 0.0), (1990.0, 0.0), (2000.0, 0.0), (2020.0, 0.0), (2100.0, 0.0)]))
    add_equation!(eqs, CEM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
    add_equation!(eqs, CIRA ~ (1 - CRDRA)^NDRA)
    add_equation!(eqs, CO2AFL ~ FOLA * (CO2AFLH / 1000) * CO2ELY * WELY)
    add_equation!(eqs, CO2AFLH ~ 1.6 * FAM)
    add_equation!(eqs, CO2ELULUC ~ CO2RFC - CO2AFL - ECO2ARA)
    add_equation!(eqs, CO2ELY ~ IfElse.ifelse(t > 2022, 1 + CO2CEACY * (CO2CA / CO2C2022 - 1), 1))
    add_equation!(eqs, CO2RFC ~ ((OGRE + CREX) * CO2RHFC) / 1000)
    add_equation!(eqs, COFE ~ FEUS * CTF / 1000)
    add_equation!(eqs, COFO ~ AFGDP * GDP + CRA + COFE)
    add_equation!(eqs, CRA ~ (ECRA * RAA) / 1000)
    add_equation!(eqs, CRBA ~ CRUS / DCS)
    add_equation!(eqs, CRBI ~ BIUS * TCTB)
    add_equation!(eqs, CRDE ~ (TUCERM + FERM + CRBI))
    add_equation!(eqs, CREX ~ IfElse.ifelse(FOLA > 0, CRLA * CREXR, 0) * ALFL * CEM)
    add_equation!(eqs, CREXR ~ 1 / 200 + CBECLE * (PCB - 1))
    add_equation!(eqs, D(CRLA) ~ CREX - CRLO - UREX)
    add_equation!(eqs, CRLO ~ CRLA * LER)
    add_equation!(eqs, CRSU ~ ACY * CRLA)
    add_equation!(eqs, CRUS ~ CRSU * (1 + CWR))
    add_equation!(eqs, CRUSP ~ CRUS / POP)
    add_equation!(eqs, CSQCA ~ ROCSQCA * SQICA)
    add_equation!(eqs, CSRA ~ CYRA * CRLA * FRA)
    add_equation!(eqs, CWR ~ ramp(t, GCWR / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, DCS ~ CRDE)
    add_equation!(eqs, DCSCA ~ DCS - CSRA)
    add_equation!(eqs, DCYCA ~ DCSCA / (CRLA * (1 - FRA)))
    add_equation!(eqs, DRM ~ ((POP * DRMP) / 1000) * (1 - FNRM))
    add_equation!(eqs, DRMP ~ TURMP)
    add_equation!(eqs, ECFT ~ CRA - FCR)
    add_equation!(eqs, ECFTSGDP ~ ECFT / GDP)
    add_equation!(eqs, ECO2ARA ~ RAA * CO2ARA / 1000)
    add_equation!(eqs, ECRA ~ ECRA22 * CIRA)
    add_equation!(eqs, FAM ~ IfElse.ifelse(t > 2022, 1 + SSP2LMA * ramp(t, (MFAM - 1) / 78, 2022, 2100), 1))
    add_equation!(eqs, FCR ~ (AFSRA / 1000) * RAA * CTF)
    add_equation!(eqs, FEER ~ 1 + FUELER * (FUCA / SFU - 1))
    add_equation!(eqs, FERM ~ RMF * KCKRM)
    add_equation!(eqs, FEUS ~ CRLA * (1 - FRA) * FUCA / 1000)
    add_equation!(eqs, FFI ~ FOFO / FF80)
    add_equation!(eqs, FFLR ~ max(0, FOLA / FOLA1980))
    add_equation!(eqs, FFLREOGRR ~ max(1, 1 + FFLREOGRRM * (FFLR - TFFLR)))
    add_equation!(eqs, FNRM ~ ramp(t, GFNRM / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, FOFO ~ CRLA * FEUS)
    add_equation!(eqs, D(FOLA) ~ NFL - CREX)
    add_equation!(eqs, FPI ~ exp(ROCFP * (t - 1980)))
    add_equation!(eqs, FRA ~ ramp(t, GFRA / IPP, 2022, 2020 + IPP))
    add_equation!(eqs, FSPI ~ exp(ROCFSP * (t - 1980)) * IfElse.ifelse(t > 2022, exp(EROCFSP * (t - 2022)), 1))
    add_equation!(eqs, FUCA ~ TFUCA / FPI)
    add_equation!(eqs, FUP ~ (FEUS / POP) * 1000)
    add_equation!(eqs, GLY ~ GLY80 + 0 * CO2CA - 0 * OW)
    add_equation!(eqs, GLY80 ~ 14 * CO2ELY * WELY)
    add_equation!(eqs, D(GRLA) ~ NGL)
    add_equation!(eqs, IUL ~ POP * ULP)
    add_equation!(eqs, LERM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
    add_equation!(eqs, LER ~ LER80 * FEER * LERM)
    add_equation!(eqs, LFL ~ OGRE + CREX)
    add_equation!(eqs, LOCR ~ CRLO + UREX)
    add_equation!(eqs, NDRA ~ log((RAA + EGB22) / EGB22) / 0.693)
    add_equation!(eqs, NFL ~ OGRE * (1 - FCG))
    add_equation!(eqs, NGL ~ OGRE * FCG)
    add_equation!(eqs, D(OGFA) ~ -NFL - NGL)
    add_equation!(eqs, OGRE ~ OGFA * OGRR * OGRRM)
    add_equation!(eqs, OGRR ~ OGRR80 * FFLREOGRR)
    add_equation!(eqs, OGRRM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
    add_equation!(eqs, PCB ~ CRBA / (1 + DRC))
    add_equation!(eqs, PRMGL ~ GRLA * GLY / 1000)
    add_equation!(eqs, RAA ~ CRLA * FRA)
    add_equation!(eqs, RMF ~ DRM - RMGL)
    add_equation!(eqs, RMGL ~ min(DRM, PRMGL))
    add_equation!(eqs, RMSP ~ (RMGL + RMF * min(1, CRBA)) * 1000 / POP)
    add_equation!(eqs, ROCSQCA ~ 0 + FUESQ * (FUCA / SFU - 1))
    add_equation!(eqs, D(SQICA) ~ CSQCA)
    add_equation!(eqs, TFA ~ OGFA + FOLA)
    add_equation!(eqs, TFUCA ~ withlookup(DCYCA, [(1.0, 0.0), (2.0, 40.0), (2.5, 50.0), (3.0, 60.0), (3.5, 70.0), (4.5, 100.0), (6.5, 200.0), (10.0, 600.0)]))
    add_equation!(eqs, TUC ~ TUCP * POP / 1000)
    add_equation!(eqs, TUCERM ~ (TUC - TUFRM) / FSPI)
    add_equation!(eqs, TUCERMP ~ TUCERM * 1000 / POP)
    add_equation!(eqs, TUCP ~ withlookup(GDPP, [(0.0, 400.0), (6.1, 680.0), (8.7, 780.0), (13.9, 950.0), (20.0, 1050.0), (30.0, 1150.0), (40.0, 1250.0), (60.0, 1350.0), (100.0, 1550.0)]))
    add_equation!(eqs, TUFRM ~ (((TURMP / 1000) * POP) - RMGL) * KCKRM)
    add_equation!(eqs, TURMP ~ withlookup(GDPP, [(0.0, 0.0), (6.1, 6.0), (8.8, 8.5), (14.0, 13.0), (30.0, 27.0), (40.0, 32.0), (50.0, 33.0), (100.0, 25.0)]))
    add_equation!(eqs, UREX ~ max(0, (IUL - URLA) / UDT))
    add_equation!(eqs, D(URLA) ~ UREX)
    add_equation!(eqs, WELY ~ IfElse.ifelse(t > 2022, 1 + OWEACY * (OW / OW2022 - 1), 1))

    # inventory
    add_equation!(eqs, CDDI ~ ROCDDI * DELDI)
    add_equation!(eqs, CPI ~ PRIN * IR)
    add_equation!(eqs, DEL ~ ((EPP / PPU) / (DELDI / DDI1980)) * IfElse.ifelse(t > 1984, PNIS, 1))
    add_equation!(eqs, D(DELDI) ~ CDDI)
    add_equation!(eqs, DSWI ~ 1 + INVEOSWI * (PRI / DRI - 1))
    smooth!(eqs, EPP, TPP, DAT)
    add_equation!(eqs, GDP ~ OUTP * PPU)
    add_equation!(eqs, IC ~ INV / RS)
    add_equation!(eqs, D(INV) ~ OUTP - DEL)
    add_equation!(eqs, IR ~ INVEOIN * (PRI / MRIWI - 1))
    add_equation!(eqs, NI ~ SA)
    add_equation!(eqs, OUTP ~ ORO * SSWI / SWI1980)
    add_equation!(eqs, D(PRIN) ~ CPI)
    add_equation!(eqs, PNIS ~ 1)
    smooth!(eqs, PRI, (IC / DIC), ICPT)
    add_equation!(eqs, ROCDDI ~ 0 + INVEODDI * (PRI / SRI - 1))
    smooth!(eqs, RS, DEL, SAT)
    add_equation!(eqs, SA ~ DEL * PPU)
    smooth!(eqs, SSWI, DSWI, TAS)

    # labourmarket
    add_equation!(eqs, AGIW ~ WARA * AHW)
    add_equation!(eqs, AHW ~ NHW / PFTJ)
    add_equation!(eqs, AHW1980 ~ NHW1980 / PFTJ80)
    add_equation!(eqs, AVWO ~ WAP * LPR)
    add_equation!(eqs, CECLR ~ ROCECLR * ECLR)
    add_equation!(eqs, CHWO ~ (OPWO - WF) / HFD)
    add_equation!(eqs, CWSO ~ WSO * ROCWSO)
    add_equation!(eqs, CWRA ~ WARA * ROCWSO)
    add_equation!(eqs, D(ECLR) ~ CECLR)
    add_equation!(eqs, GDPPEROCCLR ~ max(0, 1 + GDPPEROCCLRM * (GDPP / GDPP1980 - 1)))
    add_equation!(eqs, ENLPR2022 ~ ramp(t, GENLPR / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, HFD ~ TYLD / 3)
    add_equation!(eqs, HWMGDPP ~ 1 + TENHW * (GDPP / GDPP1980 - 1))
    add_equation!(eqs, IWEOCLR ~ 1 + WSOECLR * (WSO / WSO1980 - 1))
    add_equation!(eqs, ILPR ~ NLPR - PSW)
    add_equation!(eqs, LAPR ~ (OUTP * PRUN) / LAUS)
    add_equation!(eqs, LAUS ~ WF * AHW)
    add_equation!(eqs, LAUS80 ~ WF1980 * AHW1980)
    smooth!(eqs, LPR, ILPR, TELLM)
    add_equation!(eqs, LTEWSO ~ WSO * RWER)
    smooth!(eqs, NHW, NHW1980 * HWMGDPP, TAHW)
    add_equation!(eqs, NLPR ~ NLPR80 * (1 + WSOELPR * (WSO / WSO1980 - 1)) + ENLPR2022)
    add_equation!(eqs, OCLR ~ ECLR * WEOCLR)
    add_equation!(eqs, OPWO ~ (CAP / OCLR) * PFTJ)
    add_equation!(eqs, PART ~ LPR * (1 - PURA))
    add_equation!(eqs, PSW ~ AUR * (1 + PUELPR * (PURA / AUR - 1)))
    smooth!(eqs, PURA, UR, UPT)
    add_equation!(eqs, ROCECLR ~ ROCECLR80 * GDPPEROCCLR)
    add_equation!(eqs, ROCWSO ~ withlookup(PURA / AUR, [(0.0, 0.06), (0.5, 0.02), (1.0, 0.0), (1.5, -0.007), (2.0, -0.01)]))
    add_equation!(eqs, TCT ~ TYLD / 3)
    add_equation!(eqs, UNEM ~ max(0, AVWO - WF))
    add_equation!(eqs, UPT ~ TYLD / 3)
    add_equation!(eqs, UR ~ UNEM / AVWO)
    add_equation!(eqs, WAP ~ A20PA)
    add_equation!(eqs, D(WARA) ~ CWRA - WRE)
    add_equation!(eqs, WASH ~ WARA / LAPR)
    smooth!(eqs, WEOCLR, IWEOCLR, TCT)
    add_equation!(eqs, D(WF) ~ CHWO)
    add_equation!(eqs, WRE ~ WARA * WRER)
    add_equation!(eqs, WRER ~ IR * (1 - FIC))
    add_equation!(eqs, D(WSO) ~ CWSO - LTEWSO)

    # other
    add_equation!(eqs, CFETA ~ COFO + CE)
    add_equation!(eqs, CTA ~ CFETA)
    add_equation!(eqs, FB15 ~ 1 - (1 / (1 + exp(-LK * (GDPP - 14)))))
    add_equation!(eqs, IEL ~ 1 + INELOK * (INEQ / 0.5 - 1))
    add_equation!(eqs, LK ~ NK * IEL)
    add_equation!(eqs, PB15 ~ POP * FB15)
    add_equation!(eqs, RGGDPP ~ ((GDPP - PGDPP) / PGDPP) / TEGR)
    smooth!(eqs, PGDPP, GDPP, TEGR)

    # output
    add_equation!(eqs, AVCA ~ TS + FCI)
    add_equation!(eqs, CAPIS ~ CUCPIS / CTPIS)
    add_equation!(eqs, CAPUS ~ CUCPUS / CTPUS)
    add_equation!(eqs, CIPIS ~ max((INCPIS + OBSGIPIS * GDP) / COCA, 0))
    add_equation!(eqs, CIPUS ~ max((GIPC + OBSGIPUS * GDP) / COCA, 0))
    add_equation!(eqs, CDPIS ~ CPIS / LCPIS)
    add_equation!(eqs, CDPUS ~ CPUS / LCPUS)
    add_equation!(eqs, COCA ~ CC1980 * OWECC)
    add_equation!(eqs, D(CPIS) ~ CAPIS - CDPIS)
    add_equation!(eqs, D(CPUS) ~ CAPUS - CDPUS)
    add_equation!(eqs, CRR ~ CAPIS / CPIS)
    add_equation!(eqs, D(CUCPIS) ~ CIPIS - CAPIS)
    add_equation!(eqs, CUCPIS1980 ~ (CAPPIS1980 / LCPIS1980) * CTPIS * EMCUC)
    add_equation!(eqs, D(CUCPUS) ~ CIPUS - CAPUS)
    add_equation!(eqs, CUCPUS1980 ~ (CAPPUS1980 / LCPUS1980) * CTPUS * EMCUC)
    add_equation!(eqs, ECR ~ (ITFP - ETFP) * CRR)
    add_equation!(eqs, D(ETFP) ~ ECR)
    add_equation!(eqs, INCPIS ~ AVCA * FACNC)
    add_equation!(eqs, LCPUS ~ LCPUS1980)
    add_equation!(eqs, LCPUS1980 ~ 15 * OWELC)
    add_equation!(eqs, OBSGIPIS ~ IfElse.ifelse(t > 2022, USPIS2022, 0))
    add_equation!(eqs, OBSGIPUS ~ IfElse.ifelse(t > 2022, 0.01 + USPUS2022, 0.01))
    add_equation!(eqs, OWECC ~ IfElse.ifelse(t > 2022, 1 + OWECCM * (OW / OW2022 - 1), 1))
    add_equation!(eqs, OWELC ~ IfElse.ifelse(t > 2022, 1 + OWELCM * (OW / OW2022 - 1), 1))
    add_equation!(eqs, CAP ~ CPIS + CPUS)
    add_equation!(eqs, CBCEFCA ~ 1 + CBCEFRA * (CBC / CBC1980 - 1))
    add_equation!(eqs, EDE ~ TPP / OOV)
    add_equation!(eqs, EDEFCA ~ 1 + EDEFRA * (PEDE / ED1980 - 1))
    add_equation!(eqs, EDELC ~ 1 + EDELCM * (PEDE / ED1980 - 1))
    smooth!(eqs, FACNC, FRA1980 * FRACAMGDPPL * (WSOEFCA + CBCEFCA + EDEFCA) / 3, IPT)
    add_equation!(eqs, FRACAMGDPPL ~ max(FRACAM, 1 + GDPPEFRACA * (GDPP / GDPP1980 - 1)))
    add_equation!(eqs, ISGDP ~ (INCPIS + GIPC) / GDP)
    add_equation!(eqs, LCPIS ~ (LCPIS1980 * OWELC) / EDELC)
    smooth!(eqs, OGR, (ORO - OLY) / OLY, 1)
    smooth!(eqs, OLY, ORO, 1)
    add_equation!(eqs, OOV ~ ORO * PRUN)
    add_equation!(eqs, ORO ~ OO1980 * ((CPIS + CPUS) / (CAPPIS1980 + CAPPUS1980))^KAPPA * (LAUS / LAUS1980)^LAMBDA * (ETFP))
    smooth!(eqs, PEDE, EDE, TOED)
    add_equation!(eqs, WSOEFCA ~ 1 + WSOEFRA * (WASH / WSO1980 - 1))

    # population
    add_equation!(eqs, D(A0020) ~ BIRTHS - PASS20)
    add_equation!(eqs, D(A2040) ~ PASS20 - PASS40)
    add_equation!(eqs, A20PA ~ A2040 + A4060 + A60PL - OP)
    add_equation!(eqs, D(A4060) ~ PASS40 - PASS60)
    add_equation!(eqs, D(A60PL) ~ PASS60 - DEATHS)
    add_equation!(eqs, BIRTHR ~ BIRTHS / POP)
    add_equation!(eqs, BIRTHS ~ A2040 * FW * (OF / FP))
    add_equation!(eqs, CEFR ~ CMFR * EFR)
    add_equation!(eqs, DEATHR ~ DEATHS / POP)
    delay_n!(eqs, PASS60, RT_DEATHS, LV_DEATHS, LE60, ORDER)
    add_equation!(eqs, DEATHS ~ RT_DEATHS[ORDER])
    add_equation!(eqs, DNC ~ ((DNCM + (DNC80 - DNCM) * exp(-DNCG * (EGDPP - GDPP1980))) * (1 + DNCA * (EGDPP - GDPP1980))) * (1 - EFR) * FM)
    add_equation!(eqs, DR ~ (A0020 + A60PL) / (A2040 + A4060))
    add_equation!(eqs, EFR ~ ramp(t, GEFR / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, EPA ~ ramp(t, (GEPA - EPA22) / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, D(EGDPP) ~ (GDPP - EGDPP) / TAHI)
    add_equation!(eqs, FM ~ IfElse.ifelse(SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MFM - 1) / 78, 2022, 2100), 1), 1))
    add_equation!(eqs, GDPP ~ GDP / POP)
    add_equation!(eqs, LE ~ ((LEMAX - (LEMAX - LE1980) * exp(-LEG * (EGDPP - GDPP1980))) * (1 + LEA * (EGDPP - GDPP1980))) * WELE * LEM)
    add_equation!(eqs, LE60 ~ LE - 60)
    add_equation!(eqs, LEM ~ IfElse.ifelse(SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MLEM - 1) / 78, 2022, 2100), 1), 1))
    add_equation!(eqs, OF ~ DNC * FADFS)
    add_equation!(eqs, OP ~ A60PL * (LE - PA) / (LE - 60))
    add_equation!(eqs, PA ~ IfElse.ifelse(LE < LE1980, PA1980, PA1980 + LEEPA * (LE + EPA - LE1980)))
    delay_n!(eqs, BIRTHS, RT_PASS20, LV_PASS20, 20, ORDER)
    add_equation!(eqs, PASS20 ~ RT_PASS20[ORDER])
    delay_n!(eqs, PASS20, RT_PASS40, LV_PASS40, 20, ORDER)
    add_equation!(eqs, PASS40 ~ RT_PASS40[ORDER])
    delay_n!(eqs, PASS40, RT_PASS60, LV_PASS60, 20, ORDER)
    add_equation!(eqs, PASS60 ~ RT_PASS60[ORDER])
    add_equation!(eqs, PGR ~ BIRTHR - DEATHR)
    add_equation!(eqs, POP ~ A0020 + A2040 + A4060 + A60PL)
    add_equation!(eqs, PW ~ OP / A20PA)
    add_equation!(eqs, WELE ~ IfElse.ifelse(t > 2022, max(0, 1 + OWELE * (OW / OW2022 - 1)), 1))

    # public
    add_equation!(eqs, CTFP ~ RTA * TFPEE5TA)
    add_equation!(eqs, DRTA ~ (DROTA1980 + IfElse.ifelse(t > 2022, EDROTA2022, 0) * (1 + SCROTA * ((SC / SC1980) - 1))))
    add_equation!(eqs, ECTAF2022 ~ max(0, CTA - CTA2022))
    add_equation!(eqs, ECTAGDP ~ ECTAF2022 / GDP)
    add_equation!(eqs, GSSGDP ~ GS / GDP)
    add_equation!(eqs, IPR ~ CPUS / GPU)
    add_equation!(eqs, IROTA ~ IfElse.ifelse(t > 2022, max(0, MIROTA2022 * (1 - 1 * (GDPP / GDPTL - 1))), 0))
    add_equation!(eqs, ITFP ~ TFPEE5TA * OWTFP)
    add_equation!(eqs, OWTFP ~ IfElse.ifelse(t > 2022, 1 + OWETFP * (OW / OW2022 - 1), 1))
    add_equation!(eqs, PLUA ~ ECTAGDP * FUATA)
    add_equation!(eqs, PPP ~ max(0, 1 + IPRVPSS * log(IPR / IPR1980)))
    add_equation!(eqs, PSEP ~ VPSS / POP)
    add_equation!(eqs, PSP ~ GS / POP)
    add_equation!(eqs, D(TFPEE5TA) ~ CTFP)
    add_equation!(eqs, RROTAI ~ min(1, 1 + IIEEROTA * (INEQI / 1 - 1)))
    add_equation!(eqs, RTA ~ (DRTA + 0) * RROTAI + IROTA)
    smooth!(eqs, RTFPUA, PLUA, IPT + CTPIS)
    add_equation!(eqs, SC ~ VPSS / GDP)
    add_equation!(eqs, TFPIE5TA ~ TFPEE5TA * (1 - RTFPUA))
    add_equation!(eqs, VPSS ~ GPU * PPP)

    # wellbeing
    add_equation!(eqs, AWBDI ~ exp(DRDI + log(WDI / TDI)))
    add_equation!(eqs, AWBIN ~ 1 + IEAWBIF * (INEQ / TI - 1))
    add_equation!(eqs, AWBGW ~ max(MWBGW, min(1, 1 + GWEAWBGWF * (PWA / TW - 1))))
    add_equation!(eqs, AWBI ~ (0.5 * AWBDI + 0.5 * AWBPS) * AWBIN * AWBGW * AWBP)
    add_equation!(eqs, AWBP ~ (1 + PREAWBF * (ORP - TPR)) * WBEP)
    add_equation!(eqs, AWBPS ~ exp(DRPS + log(PSP / TPS)))
    add_equation!(eqs, IEST ~ WorldDynamics.interpolate(INEQ / AI, tables[:IEST], ranges[:IEST]))
    add_equation!(eqs, IPP ~ IfElse.ifelse(EIPF > 0, EIP, RD))
    add_equation!(eqs, IRD ~ NRD * STRERD * STEERD)
    add_equation!(eqs, IST ~ PSESTR * IEST)
    smooth!(eqs, ORP, ((AWBI - PAWBI) / AWBI) / AWBPD, AWBPD)
    smooth!(eqs, PAWBI, AWBI, AWBPD)
    add_equation!(eqs, PSSGDP ~ PSP / GDPP)
    add_equation!(eqs, PSESTR ~ WorldDynamics.interpolate(PSSGDP / SPS, tables[:PSESTR], ranges[:PSESTR]))
    smooth!(eqs, RD, IRD, TCRD)
    add_equation!(eqs, STE ~ 1 + PESTF * (ORP - AP))
    add_equation!(eqs, STEERD ~ 1 + STEERDF * (STE / STE1980 - 1))
    add_equation!(eqs, STRERD ~ 1 + STRERDF * (SOTR / STR1980 - 1))
    smooth!(eqs, SOTR, IST, TEST)
    add_equation!(eqs, WBEP ~ 1 + PAEAWBF * (LPR / THPA - 1))

    return ODESystem(eqs; name=name)
end

end;

# ╔═╡ d60dae6d-d0b7-400f-ae9e-e835100f98f7
begin
	
global p_name = ["ERDN2OKF2022", "ERDCH4KC2022", "DACCO22100", "FGDC2022", "EETF2022", "EGTRF2022", "EPTF2022", "ETGBW", "GEIC", "FETPO", "GFCO2SCCS", "EROCEPA2022", "GFNE", "GREF", "GCWR", "GFNRM", "GFRA", "EROCFSP", "USPIS2022", "USPUS2022", "GEFR", "MIROTA2022"]
global p_desc = ["Extra rate of decline in N2O per kg fertilizer from 2022", "Extra rate of decline in CH4 per kg crop after 2022 1/y", "Direct air capture of CO2 in 2100 GtCO2/y", "Fraction of government debt cancelled in 2022 1/y", "Extra empowerment tax from 2022 (share of NI)", "Extra general tax rate from 2022", "Extra pension tax from 2022 (share of NI)", "Extra transfer of government budget to workers", "Goal for extra income from commons (share of NI)", "Fraction of extra taxes paid by owners", "Goal for fraction of CO2-sources with CCS", "Extra ROC in energy productivity after 2022 1/y", "Goal for fraction new electrification", "Goal for renewable el fraction", "Goal for crop waste reduction", "Goal for fraction new red meat", "Goal for fraction regenerative agriculture", "Extra ROC in food sector productivity from 2022 1/y", "Unconventional stimulus in PIS from 2022 (share of GDP)", "Unconventional stimulus in PUS from 2022 (share of GDP)", "Goal for extra fertility reduction", "Max imported ROTA from 2022 1/y"]
global p_map = [3, 5, 6, 7, 2, 1, 12, 8, 10, 4, 15, 21, 9, 11, 13, 16, 17, 14, 22, 19, 20]
global tltl_alpha = fill(0.0, length(p_name))
global gl_alpha = fill(1.0, length(p_name))
global v_name = ["Average wellbeing", "GDP per person", "Inequality", "Observed warming", "Population", "Social tension"]

function _variables(e4a)
    variables = [
        (e4a.POP, 0, 10000, "Population"),
        (e4a.AWBI, 0, 5.0, "Average wellbeing"),
        (e4a.GDPP, 0, 65, "GDP per person"),
        (e4a.STE, 0, 2, "Social tension"),
        (e4a.INEQ, 0, 1.6, "Inequality"),
        (e4a.OW, 0, 4, "Global warming"),]
    return variables
end

function default_tltl_pars(x)
    y = copy(x)
    # Poverty turnaround
    y[4] = 0.0
    y[19] = 0.0
    y[20] = 0.0
    y[22] = 0.0
    # Inequality turnaround
    y[6] = 0.0
    y[10] = 0.5
    y[8] = 0.0
    y[9] = 0.0
    # Empowerment turnaround
    y[21] = 0.0
    y[5] = 0.0
    y[7] = 0.00
    # Food turnaround
    y[18] = 0.0
    y[15] = 0.05
    y[17] = 0.1
    y[16] = 0.1
    # Energy turnaround
    y[12] = 0.002
    y[13] = 0.5
    y[14] = 0.5
    y[11] = 0.2
    y[3] = 0.0
    # Other
    y[2] = 0.0
    y[1] = 0.0
    return y
end

function default_gl_pars(x)
    y = copy(x)
    # Poverty turnaround
    y[4] = 0.1
    y[19] = 0.01
    y[20] = 0.01
    y[22] = 0.005
    # Inequality turnaround
    y[6] = 0.01
    y[10] = 0.8
    y[8] = 0.2
    y[9] = 0.02
    # Empowerment turnaround
    y[21] = 0.2
    y[5] = 0.02
    y[7] = 0.02
    # Food turnaround
    y[18] = 0.0
    y[15] = 0.2
    y[17] = 0.5
    y[16] = 0.5
    # Energy turnaround
    y[12] = 0.004
    y[13] = 1.0
    y[14] = 1.0
    y[11] = 0.9
    y[3] = 8.0
    # Other
    y[2] = 0.01
    y[1] = 0.01
    return y
end

function e4a_sys_prob()
    @named e4a = earth4all()
    e4a_sys = structural_simplify(e4a)
    prob = ODEProblem(e4a_sys, [], (1980, 2100))
    return e4a, e4a_sys, prob
end

function compute_gl_sol(prob)
    prob = remake(prob, p=default_gl_pars(prob.p))
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function compute_tltl_sol(prob)
    prob = remake(prob, p=default_tltl_pars(prob.p))
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function plot_variables(solution, xrange, variables::Vector{<:NTuple{4,Any}}; title="", showaxis=true, showlegend=true, linetype="lines", colored=true, save=false)
    numvars = length(variables)

    @assert 1 ≤ numvars
    @assert (1 == length(xrange)) || (3 == length(xrange))
    @assert 4 == length(variables[1])


    colors = colored ? ColorSchemes.tab10.colors : fill(RGB(0.2, 0.2, 0.2), numvars)

    x_offset = 0.05
    x_domain = showaxis ? x_offset * numvars - 0.04 : 0.0


    traces = GenericTrace[]

    (xvalue, xmin, xmax) = (3 == length(xrange)) ? xrange : (xrange[1], -Inf, Inf)
    (var, varmin, varmax, varname) = variables[1]

    layout = Dict([
        ("title", attr(text=title, x=0.5)),
        ("showlegend", showlegend),
        ("plot_bgcolor", "#EEE"),
        ("xaxis", attr(
            domain=[x_domain + 0.02, 1.0],
            position=0.0,
            range=[xmin, xmax])),
        ("yaxis", attr(
            color=colors[1],
            visible=showaxis,
            name="",
            position=0.0,
            showgrid=false,
            range=[varmin, varmax],
            domain=[0.05, 1.0],
        ))
    ])

    push!(traces, scatter(
        x=solution[xvalue],
        y=solution[var],
        marker_color=colors[1],
        name=varname,
        mode=linetype, yaxis="y1"),
    )


    for i ∈ 2:numvars
        (var, varmin, varmax, varname) = variables[i]

        layout[string("yaxis", i)] = attr(
            color=colors[i],
            overlaying="y",
            visible=showaxis,
            name="",
            position=(i - 1) * x_offset,
            showgrid=false,
            range=[varmin, varmax],
        )

        push!(traces, scatter(
            x=solution[xvalue],
            y=solution[var],
            marker_color=colors[i],
            name=varname,
            mode=linetype,
            yaxis=string("y", i)),
        )
    end

    p = plot(traces, Layout(layout))
    save && savefig(p, "./" * title * ".svg")

    return p
end

function plot_sol_var(e4a, _sol, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_variables(_sol, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
end

function plot_two_sol_variables(scen1, sol1, scen2, sol2, xrange, variables::Vector{<:NTuple{4,Any}}; title="", showaxis=true, showlegend=true, linetype="lines", colored=true, save=false)
    numvars = length(variables)

    @assert 1 ≤ numvars
    @assert (1 == length(xrange)) || (3 == length(xrange))
    @assert 4 == length(variables[1])


    colors = colored ? ColorSchemes.tab10.colors : fill(RGB(0.2, 0.2, 0.2), numvars)

    x_offset = 0.05
    x_domain = showaxis ? x_offset * numvars - 0.04 : 0.0


    traces = GenericTrace[]

    (xvalue, xmin, xmax) = (3 == length(xrange)) ? xrange : (xrange[1], -Inf, Inf)
    (var, varmin, varmax, varname) = variables[1]

    layout = Dict([
        ("title", attr(text=title, x=0.5)),
        ("showlegend", showlegend),
        ("plot_bgcolor", "#EEE"),
        ("xaxis", attr(
            domain=[x_domain + 0.02, 1.0],
            position=0.0,
            range=[xmin, xmax])),
        ("yaxis", attr(
            color=colors[1],
            visible=showaxis,
            name="",
            position=0.0,
            showgrid=false,
            range=[varmin, varmax],
            domain=[0.05, 1.0],
        ))
    ])

    push!(traces, scatter(
        x=sol1[xvalue],
        y=sol1[var],
        marker_color=colors[1],
        name=varname * " (" * scen1 * ")",
        mode=linetype, yaxis="y1"),
    )

    push!(traces, scatter(
        x=sol2[xvalue],
        y=sol2[var],
        marker_color=colors[1],
        name=varname * " (" * scen2 * ")",
        mode=linetype, yaxis="y1",
        line=attr(dash="dash"))
    )

    for i ∈ 2:numvars
        (var, varmin, varmax, varname) = variables[i]

        layout[string("yaxis", i)] = attr(
            color=colors[i],
            overlaying="y",
            visible=showaxis,
            name="",
            position=(i - 1) * x_offset,
            showgrid=false,
            range=[varmin, varmax],
        )

        push!(traces, scatter(
            x=sol1[xvalue],
            y=sol1[var],
            marker_color=colors[i],
            name=varname * " (" * scen1 * ")",
            mode=linetype,
            yaxis=string("y", i)),
        )

        push!(traces, scatter(
            x=sol2[xvalue],
            y=sol2[var],
            marker_color=colors[i],
            name=varname * " (" * scen2 * ")",
            mode=linetype,
            yaxis=string("y", i),
            line=attr(dash="dash"))
        )
    end

    p = plot(traces, Layout(layout))
    save && savefig(p, "./" * title * ".svg")

    return p
end

function plot_two_sol_var(e4a, _scen1, _sol1, _scen2, _sol2, var_ind, _title, sa, sl; kwargs...)
    @variables t
    return plot_two_sol_variables(_scen1, _sol1, _scen2, _sol2, (t, 1980, 2100), getindex(_variables(e4a), var_ind); title=_title, showaxis=sa, showlegend=sl, kwargs...)
end

function compute_alpha_sol(prob, fixed_p, α, p)
    tltl_p = default_tltl_pars(prob.p)
    gl_p = default_gl_pars(prob.p)
    δ = gl_p - tltl_p
    for i in 1:lastindex(p)
        tltl_p[p[i]] = tltl_p[p[i]] + δ[p[i]] * α[i]
    end
    for i in 1:lastindex(fixed_p)
        tltl_p[fixed_p[i][1]] = fixed_p[i][2]
    end
    prob = remake(prob, p=tltl_p)
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function compute_gl18_sol(prob)
    a = fill(1.0, 22)
    a[1] = 0.0
    a[12] = 0.0
    a[15] = 0.0
    return compute_alpha_sol(prob, [], a, collect(1:22))
end

	
function plot_two_sols(scen1, sol1, scen2, sol2, vars, plot_title)
    x = range(1, 7681, length=7681)
    traces = []
    lbl = fill("NA", 1, 2*length(vars))
    for v in 1:lastindex(vars)
		desc = getdescription(vars[v])
		push!(traces, sol1[vars[v]])
		push!(traces, sol2[vars[v]])
		lbl[1,2*v-1] = getdescription(vars[v]) * "-" * scen1
		lbl[1,2*v] = getdescription(vars[v]) * "-" * scen2
    end
    return plot(x, traces, xticks = ([1:640:7681;],["1980","1990","2000","2010","2020","2030","2040","2050","2060","2070","2080","2090","2100"]), label=lbl, title=plot_title, legend=:outertopright)
end

function compute_sol(prob, pluto_p)
	perm_p = prob.p
	for i in 1:lastindex(pluto_p)
		perm_p[p_map[i]] = pluto_p[i]
	end
    prob = remake(prob, p=perm_p)
    sol = DifferentialEquations.solve(prob, Euler(); dt=0.015625, dtmax=0.015625)
    return sol
end

function plot_three_sols(scen1, sol1, scen2, sol2, scen3, sol3, vars, plot_title)
    x = range(1, 7681, length=7681)
    traces = []
    lbl = fill("NA", 1, 3*length(vars))
    for v in 1:lastindex(vars)
		desc = getdescription(vars[v])
		push!(traces, sol1[vars[v]])
		push!(traces, sol2[vars[v]])
		push!(traces, sol3[vars[v]])
		lbl[1,3*v-2] = getdescription(vars[v]) * "-" * scen1
		lbl[1,3*v-1] = getdescription(vars[v]) * "-" * scen2
		lbl[1,3*v] = getdescription(vars[v]) * "-" * scen3
    end
    return plot(x, traces, xticks = ([1:640:7681;],["1980","1990","2000","2010","2020","2030","2040","2050","2060","2070","2080","2090","2100"]), label=lbl, title=plot_title, legend=:outertopright)
end


end;

# ╔═╡ 2a90243e-2457-49fa-a8b9-7e28d6b4d4ca
md"""
We can now compute the seven scenarios analyzed in the paper.
"""

# ╔═╡ c55f68b4-fa85-4e54-874f-811810282fdf
begin
	e4a, _, prob = e4a_sys_prob()
	tltl = compute_tltl_sol(prob)
	gl = compute_gl_sol(prob)
	gl18 = compute_gl18_sol(prob);
	gl6 = compute_alpha_sol(prob, [], [1.0,1.0,1.0,1.0,1.0,1.0], [3,9,11,19,20,22])
	dgl6 = compute_alpha_sol(prob, [], [2.0,2.0,2.0,2.0,2.0,2.0], [3,9,11,19,20,22])
	gl4 = compute_alpha_sol(prob, [], [0.25,2.0,0.5,2.0,2.0,2.0], [3,9,11,19,20,22])
	a = fill(1.0, 22)
	a[3] = 0.25
	a[11] = 0.5
	mgl = compute_alpha_sol(prob, [], a, collect(1:22))
	println("All scenarios have been computed")
end

# ╔═╡ b2693f34-b19a-42b2-98a5-d112b53bc48c
md"""
Below we can select two scenarios (among the computed seven ones) and one output indicator, in order to compare the evolution of the indicator in the two scenarios. 
"""

# ╔═╡ 03833b31-f9a3-49ea-971b-e1b83a259b0f
md"""
First scenario $(@bind scen1 Select(["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]))
Second scenario $(@bind scen2 Select(["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]))
Indicator $(@bind indicator Select([e4a.AWBI,e4a.GDPP,e4a.INEQ,e4a.OW,e4a.POP,e4a.STE]))
"""

# ╔═╡ 4e516d92-2e77-46e1-a8a4-7377a5912caf
begin
	plotly()
	scens = [tltl,gl,gl18,gl6,dgl6,gl4,mgl]
	scen_names = ["TLTL","GL","GL18","GL6","DGL6","GL4","MGL"]
	i1 = findfirst(s -> s==scen1, scen_names)
	i2 = findfirst(s -> s==scen2, scen_names)
	plot_two_sols(scen1, scens[i1], scen2, scens[i2], [indicator], scen1 * " versus " * scen2)
end


# ╔═╡ c117bbaf-11be-4c58-a9cd-843593b77729
md"""
In order to compare the evolution of the other variables of the model, we can specify the acronym of the variable in the code below prefized by `e4a`.
"""

# ╔═╡ 22d18146-71ab-4838-8c0d-81c91a06951e
plot_two_sols(scen1, scens[i1], scen2, scens[i2], [e4a.OF], scen1 * " versus " * scen2)

# ╔═╡ cca46575-3452-4aca-b80f-c793c3c3dcd4
md"""
Finally, below we can set the value of the twenty-one parameters involved in the five turnarounds and (by specifying its acronym prefixed by `e4a`) see the evolution of any variable in the TLTL scenario, in the GL scenario and in the scenario (called Pluto) produced by the chosen set of values.
"""

# ╔═╡ 1f1439a4-dc73-4265-b264-10c69ec311ba
begin
struct MySlider 
    range::AbstractRange
    default::Number
end
function Base.show(io::IO, ::MIME"text/html", slider::MySlider)
    print(io, """
		<input type="range" 
		min="$(first(slider.range))" 
		step="$(PlutoUI.step(slider.range))"
		max="$(last(slider.range))" 
		value="$(slider.default)"
		oninput="this.nextElementSibling.value=this.value">
		<output>$(slider.default)</output>""")
end
end

# ╔═╡ c3ec6930-816e-41b7-94c7-50a862216f27
begin
	md"Direct air capture of CO2 in 2100 GtCO2/y $(@bind DACCO22100 MySlider(0.0:0.5:16.0,8.0))\
	Extra empowerment tax from 2022 $(@bind EETF2022 MySlider(0.0:0.001:0.04,0.02))\
	Extra general tax rate from 2022 $(@bind EGTRF2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra pension tax from 2022 (share of NI) $(@bind EPTF2022 MySlider(0.0:0.001:0.04,0.02))\
	Extra rate of decline in CH4 per kg crop after 2022 1/y $(@bind ERDCH4KC2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra rate of decline in N2O per kg fertilizer from 2022 $(@bind ERDN2OKF2022 MySlider(0.0:0.001:0.02,0.01))\
	Extra ROC in energy productivity after 2022 1/y $(@bind EROCEPA2022 MySlider(0.002:0.0001:0.006,0.004))\
	Extra transfer of govmnt budget to workers $(@bind ETGBW MySlider(0.0:0.01:0.4,0.2))\
	Fraction of extra taxes paid by owners $(@bind FETPO MySlider(0.5:0.01:0.99,0.8))\
	Fraction of govmnt debt cancelled in 2022 1/y $(@bind FGDC2022 MySlider(0.0:0.01:0.2,0.1))\
	Goal for crop waste reduction $(@bind GCWR MySlider(0.05:0.01:0.35,0.2))\
	Goal for extra fertility reduction $(@bind GEFR MySlider(0.0:0.01:0.4,0.2))\
	Goal for extra income from commons (share of NI) $(@bind GEIC MySlider(0.0:0.001:0.04,0.02))\
	Goal for fraction of CO2-sources with CCS $(@bind GFCO2SCCS MySlider(0.2:0.01:0.99,0.9))\
	Goal for fraction new electrification $(@bind GFNE MySlider(0.5:0.01:1.0,1.0))\
	Goal for fraction new red meat $(@bind GFNRM MySlider(0.1:0.01:0.9,0.5))\
	Goal for fraction regenerative agriculture $(@bind GFRA MySlider(0.1:0.01:0.9,0.5))\
	Goal for renewable el fraction $(@bind GREF MySlider(0.5:0.01:1.0,1.0))\
	Max imported ROTA from 2022 1/y $(@bind MIROTA2022 MySlider(0.0:0.0001:0.01,0.005))\
	Unconventional stimulus in PIS from 2022 (share of GDP) $(@bind USPIS2022 MySlider(0.0:0.001:0.02,0.01))\
	Unconventional stimulus in PUS from 2022 (share of GDP) $(@bind USPUS2022 MySlider(0.0:0.001:0.02,0.01))\
	"
end


# ╔═╡ 9105a4ed-e4d4-42d0-95b4-90093cbf8c1e
sol = compute_sol(prob, [DACCO22100,EETF2022,EGTRF2022,EPTF2022,ERDCH4KC2022,ERDN2OKF2022,EROCEPA2022,ETGBW,FETPO,FGDC2022,GCWR,GEFR,GEIC,GFCO2SCCS,GFNE,GFNRM,GFRA,GREF,MIROTA2022,USPIS2022,USPUS2022]);


# ╔═╡ 9e4c3d9d-7249-421f-a7fe-b0221322324c
plot_three_sols("TLTL", tltl, "GL", gl, "Pluto", sol, [e4a.STE], "Pluto versus TLTL and GL")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
IfElse = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
WorldDynamics = "34572acc-64b2-43ec-bac8-7b2a4baad13f"

[compat]
DifferentialEquations = "~7.12.0"
IfElse = "~0.1.1"
ModelingToolkit = "~8.75.0"
Plots = "~1.40.3"
PlutoUI = "~0.7.58"
WorldDynamics = "~0.4.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "96e0e744b03f63d523beca5d080854186b61e943"

[[deps.ADTypes]]
git-tree-sha1 = "016833eb52ba2d6bea9fcb50ca295980e728ee24"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.7"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Markdown", "Test"]
git-tree-sha1 = "c0d491ef0b135fd7d63cbc6404286bc633329425"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.36"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"
    AccessorsUnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "44691067188f6bd1b2289552a23e4b7572f4528d"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.9.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "6404a564c24a994814106c374bec893195e19bac"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.8.0"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "c946c5014cf4cdbfacacb363b110e7bffba3e742"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.6.1"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "c9b163bd832e023571e86d0b90d9de92a9879088"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.6"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.Blink]]
deps = ["Base64", "Distributed", "HTTP", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Pkg", "Reexport", "Sockets", "WebIO"]
git-tree-sha1 = "bc93511973d1f949d45b0ea17878e6cb0ad484a1"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.9"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "BandedMatrices", "ConcreteStructs", "DiffEqBase", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "005b55fa2eebaa4d7bf3cfb8097807f47116175f"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.7.1"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Static"]
git-tree-sha1 = "601f7e7b3d36f18790e2caf83a882d88e9b71ff1"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.4"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "b544d62417a99d091c569b95109bc9d8c223e9e3"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.4.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "575cd02e080939a33b6df6c5853d14924c08e35b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.23.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "0f4b5d62a88d8f59003e43c25a8a90de9eb76317"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.18"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "bfae672496149b369172eae6296290a381df2bdf"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.47.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces"]
git-tree-sha1 = "4fa023dbb15b3485426bbc6c43e030c14250d664"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.149.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "Functors", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "ee954c8b9d348b7a8a6aec5f28288bf5adecd4ee"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.37.0"
weakdeps = ["OrdinaryDiffEq", "Sundials"]

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "65cbbe1450ced323b4b17228ccd96349d96795a7"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.21.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "8864b6a953eeba7890d23258aca468d90ca73fd6"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.12.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "51b4b84d33ec5e0955b55ff4b748b99ce2c3faa9"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.7"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "0bb0a6f812213ecc8fbbcf472f4a993036858971"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.5.5"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.EnzymeCore]]
git-tree-sha1 = "1bc328eec34ffd80357f84a84bb30e4374e9bd60"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.6.6"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "8e18940a5ba7f4ddb41fe2b79b6acaac50880a86"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.26.1"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "9dc913faf8552fd09b92a0d7fcc25f1d5609d795"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "a6e756a880fc419c8b41592010aebe6a5ce09136"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.8"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0a59c7d1002f3131de53dc4568a47d15a44daef7"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "2.0.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "bc0c5092d6caaea112d3c8e3b238d61563c58d5f"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fa8d8fcfa6c38a9a7aa07233e35b3d9a39ec751a"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.9"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3437ade7073682993e092ca570ad68a2aba26983"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a96d5c713e6aa28c242b0d25c1347e258d6541ab"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.3+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "af49a0851f8113fcfae2ef5027c6d49d0acec39b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.4"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

[[deps.Glob]]
git-tree-sha1 = "97285bbd5230dd766e9ef6749b80fc617126d496"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "eb8fed28f4994600e29beef49744639d985a04b2"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.16"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "896385798a8d49a255c398bd49162062e4a4c435"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.13"
weakdeps = ["Dates"]

    [deps.InverseFunctions.extensions]
    DatesExt = "Dates"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Glob", "Pkg", "PrecompileTools", "Tokenize"]
git-tree-sha1 = "e07d6fd7db543b11cd90ed764efec53f39851f09"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "1.0.54"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "ed08d89318be7d625613f3c435d1f6678fba4850"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.11.1"
weakdeps = ["FastBroadcast"]

    [deps.JumpProcesses.extensions]
    JumpProcessFastBroadcastExt = "FastBroadcast"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "07649c499349dad9f08dde4243a4c597064663e9"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.6.0"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43032da5832754f58d14a91ffbe86d5f176acda9"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.2.1+0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "8a6837ec02fe5fb3def1abc907bb802ef11a0729"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "d1f981fba6eb3ec393eede4821bca3f2b7592cd4"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.15.1"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "62edfee3211981241b57ff1cedf4d74d79519277"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.15"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays"]
git-tree-sha1 = "af45931c321aafdb96a6e0b26e81124e1b390e4e"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "1.9.0"
weakdeps = ["StaticArrays"]

    [deps.LazyArrays.extensions]
    LazyArraysStaticArraysExt = "StaticArrays"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "InteractiveUtils", "KLU", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "PrecompileTools", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SparseArrays", "Sparspak", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "775e5e5d9ace42ef8deeb236587abc69e70dc455"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "2.28.0"

    [deps.LinearSolve.extensions]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveEnzymeExt = ["Enzyme", "EnzymeCore"]
    LinearSolveFastAlmostBandedMatricesExt = ["FastAlmostBandedMatrices"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = "Pardiso"
    LinearSolveRecursiveArrayToolsExt = "RecursiveArrayTools"

    [deps.LinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "PrecompileTools", "SIMDTypes", "SLEEFPirates", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "a13f3be5d84b9c95465d743c82af0b094ef9c2e2"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.169"
weakdeps = ["ChainRulesCore", "ForwardDiff", "SpecialFunctions"]

    [deps.LoopVectorization.extensions]
    ForwardDiffExt = ["ChainRulesCore", "ForwardDiff"]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "78f6e33434939b0ac9ba1df81e6d005ee85a7396"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "2.1.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "b1f2f92feb0bc201e91c155ef575bcc7d9cc3526"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.2"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "Combinatorics", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "ForwardDiff", "FunctionWrappersWrappers", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "MacroTools", "NaNMath", "OrdinaryDiffEq", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "11f35f9619c625c18454a94a9013e7d047501fbf"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.75.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKDeepDiffsExt = "DeepDiffs"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "769c9175942d91ed9b83fa929eee4fe6a1d128ad"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.4"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "a7cefa21a2ff993bff0456bf7521f46fc077ddf1"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.19"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "2d106538aebe1c165e16d277914e10c550e9d9b7"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.2"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "MbedTLS", "Pkg", "Sockets"]
git-tree-sha1 = "0bdaa479939d2a1f85e2f93e38fbccfcb73175a5"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "1.0.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "LazyArrays", "LineSearches", "LinearAlgebra", "LinearSolve", "MaybeInplace", "PrecompileTools", "Preferences", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "TimerOutputs"]
git-tree-sha1 = "b9e12aa04c90a05d2aaded6f7c4d8b39e77751db"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "3.9.1"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveBandedMatricesExt = "BandedMatrices"
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = "NLsolve"
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSymbolicsExt = "Symbolics"
    NonlinearSolveZygoteExt = "Zygote"

    [deps.NonlinearSolve.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d9b79c4eed437421ac4285148fcadf42e0700e89"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.4"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "91079af18db922354197eeae2a17b177079e24c1"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.74.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "bf9f8fbb4e76a621fe3a54a7281c5410b77ba015"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.12"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "3bdfa4fa528ef21287ef659a89d686e8a1bcb1a9"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.3"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "09f59c6dda37c7f73efddc5bdf6f92bc940eb484"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.12"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff"]
git-tree-sha1 = "b6665214f2d0739f2d09a17474dd443b9139784a"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.20"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "4743b43e5a9c4a2ede372de7061eed81795b12e7"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "SparseArrays", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d8f131090f2e44b145084928856a561c83f43b27"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.13.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "PrecompileTools", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "8bc86c78c7d8e2a5fe559e3721c0f9c9e303b2ed"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.21"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "6aacc5eefe8415f47b3e34214c1d79d2674a0ba2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "3aac6d68c5e57449f5b9b865c9ba50ac2970c4cf"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.42"

[[deps.SciMLBase]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "d15c65e25615272e1b1c5edb1d307484c7942824"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.31.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseMakieExt = "Makie"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseZygoteExt = "Zygote"

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools", "Setfield", "SparseArrays", "StaticArraysCore"]
git-tree-sha1 = "10499f619ef6e890f3f4a38914481cc868689cd5"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.8"

[[deps.SciMLStructures]]
git-tree-sha1 = "5833c10ce83d690c124beedfe5f621b50b02ba4d"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.1.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "DiffResults", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "Reexport", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "d4c17fc60bf5f8f2be02777c4836878f27ac7b9b"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "1.7.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolvePolyesterForwardDiffExt = "PolyesterForwardDiff"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveStaticArraysExt = "StaticArrays"
    SimpleNonlinearSolveTrackerExt = "Tracker"
    SimpleNonlinearSolveZygoteExt = "Zygote"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "a616ac46c38da60ac05cecf52064d44732edd05e"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.17.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterExt = "Polyester"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Polyester = "f517fe37-dbe3-4b94-8317-1923a5111588"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "d2fdac9ff3906e27f7a618d47b676941baa6c80c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.10"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "bf074c045d3d5ffd956fa0a461da38a44685d6b2"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.3"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "a735fd5053724cf4de31c81b4e2cc429db844be5"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.0.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "97e5d0b7e5ec2e68eec6626af97c59e9f6b6c3d0"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.65.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "d6415f66f3d89c615929af907fdc6a3e17af0d8c"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.2"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "e15f5a73f0d14b9079b807a9d1dac13e4302e997"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.24.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "ba4d38faeb62de7ef47155ed321dce40a549c305"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.2+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "MacroTools", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "4b7f4c80449d8baae8857d55535033981862619c"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.15"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "89aa6b25a75418c8fffc42073b2e7dce69847394"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "669e43e90df46fcee4aa859b587da7a7948272ac"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.5.1"

[[deps.Symbolics]]
deps = ["ArrayInterface", "Bijections", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "ForwardDiff", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "PrecompileTools", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils"]
git-tree-sha1 = "280c17e091a24283a59eedfc00a02026ec984b09"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.27.1"

    [deps.Symbolics.extensions]
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxCoreExt = "LuxCore"
    SymbolicsPreallocationToolsExt = "PreallocationTools"
    SymbolicsSymPyExt = "SymPy"

    [deps.Symbolics.weakdeps]
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    LuxCore = "bb33d45b-7691-41d6-9220-0943567d0623"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "eda08f7e9818eb53661b3deb74e3159460dfbc27"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.2"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.Tokenize]]
git-tree-sha1 = "5b5a892ba7704c0977013bd0f9c30f5d962181e0"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.28"

[[deps.TranscodingStreams]]
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "7ee8ed8904e7dd5d31bb46294ef5644d9e2e44e4"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.21"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "7209df901e6ed7489fe9b7aa3e46fb788e15db85"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.65"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WorldDynamics]]
deps = ["ColorSchemes", "ColorTypes", "DifferentialEquations", "IfElse", "ModelingToolkit", "PlotlyJS"]
git-tree-sha1 = "060e2bbc95c527ce3e1ee91180937eb4f0ccd6d8"
uuid = "34572acc-64b2-43ec-bac8-7b2a4baad13f"
version = "0.4.4"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "532e22cf7be8462035d092ff21fada7527e2c488"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.6+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─b7e9f8cc-f1b1-11ee-19fa-9de56c9de796
# ╟─0792fd16-759f-4c6d-8973-be41a4510cb8
# ╠═b286eea2-1258-4cf9-a14a-1eedd94bd387
# ╟─78c67880-9ba9-410b-be26-27ea66f97d4a
# ╟─fc6c03cf-a831-4967-8446-df9b8d0928d8
# ╟─0c632ed8-0fc0-4f41-a1e8-fb2a40438857
# ╟─f2d85dc2-d2db-4953-8d49-4051bc5599e8
# ╟─a30de87a-fae9-4270-b1da-3a80de97b235
# ╟─d60dae6d-d0b7-400f-ae9e-e835100f98f7
# ╟─2a90243e-2457-49fa-a8b9-7e28d6b4d4ca
# ╠═c55f68b4-fa85-4e54-874f-811810282fdf
# ╟─b2693f34-b19a-42b2-98a5-d112b53bc48c
# ╟─03833b31-f9a3-49ea-971b-e1b83a259b0f
# ╟─4e516d92-2e77-46e1-a8a4-7377a5912caf
# ╟─c117bbaf-11be-4c58-a9cd-843593b77729
# ╟─22d18146-71ab-4838-8c0d-81c91a06951e
# ╟─cca46575-3452-4aca-b80f-c793c3c3dcd4
# ╟─1f1439a4-dc73-4265-b264-10c69ec311ba
# ╟─c3ec6930-816e-41b7-94c7-50a862216f27
# ╟─9105a4ed-e4d4-42d0-95b4-90093cbf8c1e
# ╠═9e4c3d9d-7249-421f-a7fe-b0221322324c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
