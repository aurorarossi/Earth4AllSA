_params = Dict{Symbol,Float64}(
    :INITFOLA => 1100,
    :INITNHW => 2,
    :INITWSO => 0.5,
    :INITWF => 1530,
    :GDPP1980 => 6.4,
    :EPA22 => 0,
    :LE1980 => 67,
    :PA1980 => 62,
    :INITSTE => 1.3,
    :INITSOTR => 0.6,
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
