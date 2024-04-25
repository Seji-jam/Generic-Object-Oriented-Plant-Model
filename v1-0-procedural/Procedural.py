import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import Functions

pi=3.141592653589793


# Weather and Nitrogen States
WetIdx = -1
NitIdx = -1
WatInp = 15
NitApp = 0
NitFix = 0.65

# Plant Types and Growth Factors
LegTyp = 1
C3C4Ty = 1
GrowP = -1
SlopeV = -1
FallR = 1

# Crop Fraction Values and Nitrogen Concentration
GrEff = 0.35
CarFV = 0.473  # Carbon fraction in the vegetative parts
YielGV = 0.80
FatFr = 0.02
LignFr = 0.06
OAcidF = 0.04
MinNF = 0.03
LeafNC = 0.055  # Nitrogen concentration in living leaves

# Temperature and Growth Parameters
BaseTD = 0
OptTD = 27.6
CritTD = 36
SenTD = 0.409
StartP = 0.2
EndP = 0.7
InitSP = -2
LeafWd = 0.025
RootDm = 100
CropHt = 345
MaxEH = 0.8
DrySI = 1.35
MaxES = 0.5
CarbonF = 6
NUpTk = 0.65
LeafSA = 0.0333
MinSLN = 0.5
MinRCN = 0.005
StemNC = 0.015
RootB = 0.25
MaxAJ = 48041.88
ValN = 62
ValJ = 124
WatUse = 0.7

# Seed Weight and Nitrogen Content
SeedWt = 0.21480
SeedNC = 0.04625
BulkD = 50.0
MaxHt = 0.7054
DevDurV = 34.7627
DevDurR = 23.0889
SenSc = -0.0

# Soil Moisture and Carbon Parameters
PanLS = 1
SoilTy = 23.4
MinWC = 0.05
FieldC = 0.25
MaxWC = 0.35
DropM = 1.44
DecPM = 10.0
RootPM = 0.3
BioOM = 0.66
HumusR = 0.02
TotalOC = 7193.0
BioCH = 3500.0
FrBioC = 0.03
NRatio = 1.0
AirR = 1
SoilSS = 100.0
SeedD = 5.0
TillCT = 4.0
TillCP = 1.0
MulF = 1

# CO2 and Respiration Rates
CO2Atm = 350.0
ResRF = 1.0
ResVF = 1.0
TempF = 0
CrushF = 0.5
NitRSh = 0.63
PrePN = 0.7
CarbonB = 0.75
CarbonX = 1.0
TempM = 1.5

# Timing and Plant Numbers
StartD = 110
NumPL = 60
NumSV = 53
NumEQ = NumSV

# Arrays for States and Rates
States = np.zeros(NumSV)
Rates = np.zeros(NumSV)

# Initial conditions
FruitP = 6.25 * SeedNC
LeafNMin = LeafSA * MinSLN
RootDepthIni = max(2.0, SeedD)
InitH = MaxHt / 1000.0
StressF = 0.9

# Initial conditions for the soil model
SoilTIni = 15
WatCI = FieldC * MulF
DecPI = 0
BioIni = FrBioC * TotalOC
NAppI = 2.0
NNitI = 2.0
WUpLI = 10.0 * (WatCI - MinWC) * RootDepthIni
WLoLI = 10.0 * (WatCI - MinWC) * (150.0 - RootDepthIni)
RootPMI = TotalOC - BioCH - DecPI
HumIni = BioCH - BioIni
DPNitI = 1.0 / 40.0 * DecPI
NUpLIni = (1.0 - np.exp(-0.065 * RootDepthIni)) * NAppI + RootDepthIni / 150.0 * AirR
NLoLIni = np.exp(-0.065 * RootDepthIni) * NAppI + (1.0 - RootDepthIni / 150.0) * AirR
NNUpLIni = (1.0 - np.exp(-0.065 * RootDepthIni)) * NNitI + RootDepthIni / 150.0 * NRatio
NNLoLIni = np.exp(-0.065 * RootDepthIni) * NNitI + (1.0 - RootDepthIni / 150.0) * NRatio
CarbFR = 1.0 - FruitP - FatFr - LignFr - OAcidF - MinNF
RPNitI = 1.0 / 100.0 * RootPMI
CFormO = 0.444 * CarbFR + 0.531 * FruitP + 0.774 * FatFr + 0.667 * LignFr + 0.368 * OAcidF
LeafCIni = NumPL * SeedWt * CFormO * GrEff * CrushF  # Initial amount of carbon in the living leaves

RootCIni = NumPL * SeedWt * CFormO * GrEff * (1.0 - CrushF)
YieldGO = CFormO / (1.275 * CarbFR + 1.887 * FruitP + 3.189 * FatFr + 2.231 * LignFr + 0.954 * OAcidF) * 30.0 / 12.0
NitLVI = LeafNC * LeafCIni / CarFV
LAIIni = LeafCIni / CarFV * LeafSA
NitRTI = NumPL * SeedWt * GrEff * LeafNC * CrushF / NitRSh - NitLVI
SLNBIni = NitLVI / LAIIni

DayShift = 0
CumTDU = 0
LeafCD = 0
StemSC = 0
CropSC = 0

RootCD = 0
LeafCDS = 0

StemN = 0

CropNS = 0

LeafND = 0
RootND = 0
RootVC = 0
RootVR = 0
NitREndE = 0
NitREndF = 0
DCDSR = 0
DCDTR = 0
RspMul = 0
NitDemP = 0
NSupplyPP = 0
NitFixT = 0
DCDTP = 0
TotalDAPAR = 0
TotalAPCan = 0
TotalResp = 0
TotalTranspCan = 0
TotalNUptake = 0
LitNitroT = 0
TotalNNLoL = 0
SoilFernA = 0

SLNBIni = SLNBIni
TotalNLV = NitLVI
NitLVI = NitLVI
NitRT = NitRTI
RootSC = RootCIni
LeafCIni = LeafCIni
LAIIni = LAIIni
ElevM = InitH
RootD = RootDepthIni
SoilTIni = SoilTIni
WUpLI = WUpLI
WLoLI = WLoLI
DecPI = DecPI
RootPMI = RootPMI
BioIni = BioIni
HumIni = HumIni
DPNitI = DPNitI
RPNitI = RPNitI
NUpLIni = NUpLIni
NLoLIni = NLoLIni
NNUpLIni = NNUpLIni
NNLoLIni = NNLoLIni


DELT = 1
PRDEL = 1




NLD = pd.read_excel(r'Weather_File.xlsx', sheet_name='Sheet1').values
OutTbl = np.zeros((NLD.shape[0], 8))

for I in range(0, NLD.shape[0]):

    LatPos = 51.97
    LonPos = 5.4
    ElevM = 7

    WeatherSt = NLD[I, 0]
    Year2003 = NLD[I, 1]
    DayYr = NLD[I, 2]
    DayDegT = NLD[I, 3] * 1000
    TminD = NLD[I, 4]
    TmaxD = NLD[I, 5]
    VapPres = NLD[I, 6]
    WindSp = NLD[I, 7]
    Precip = NLD[I, 8]

    CurrDay = DayYr

    print('      >>>>>     ', I, ':     DOY =', CurrDay)

    DayShift = CurrDay - StartD + 1

    AvgTemp = 0.29 * TminD + 0.71 * TmaxD
    NAvgTemp = 0.71 * TminD + 0.29 * TmaxD
    
    # Weather and Day Length Calculations
    [SunCont, SinLatD, CosLatD, DayLen, DayDegL, DSinB] = Functions.WEATHER(DayYr, LatPos, InitSP)
    
    # Soil Water Availability
    WCUpperL = (WUpLI + MinWC * 10. * RootD) / 10. / RootD # Soil water content of the upper soil layer 
    WCLowerL = min(MaxWC, (WLoLI + MinWC * 10. * (150. - RootD)) / 10. / (150. - RootD)) # Soil water content of the lower soil layer 
    DailyWSupply = Functions.INSWTICH(WetIdx, WatInp, max(0.1, WUpLI / TillCP + 0.1)) # Daily water supply for evapotranspiration 

    LodgeS = Functions.INSWTICH(FallR, 0, 0)
    
    # Leaf Area Development
    DeltaLAI = (LeafCIni - RootCIni) / CarFV * LeafSA 
    TotalLAI = LAIIni + LeafCIni / CarFV * LeafSA
    
    # Extinction Coefficient 
    KLeaf = Functions.KDF_COEFF(TotalLAI, BulkD * 3.141592654 / 180., 0.2)
    KLeafN = KLeaf * (NitLVI - MinSLN * TotalLAI)
    NBack = MinSLN * (1.0 - np.exp(-KLeaf * TotalLAI))
    KWater = KLeaf
    KNitro = 1.0 / TotalLAI * np.log((KLeafN + NBack) / (KLeafN * np.exp(-KLeaf * TotalLAI) + NBack))
    
    # Leaf Area Development Continued
    LAIAdj = np.log(1. + KNitro * max(0., NitLVI) / MinSLN) / KNitro
    LAIFin = min(LAIAdj, LAIIni)
    # SLA  = LAIFin / Leaf Weight (WLV is not defined)

    # Specific Leaf Nitrogen 
    SLNitro = NitLVI / LAIFin  # Specific leaf nitrogen content (average value in canopy)
    SLNTop = NitLVI * KNitro / (1. - np.exp(-KNitro * LAIFin)) # SLN for top leaves in canopy
    # SLNB (Specific Leaf Nitrogen Bottom) is in bottom leaves of canopy which initially is set to zero
    SLNBtmC = NitLVI * KNitro * np.exp(-KNitro * LAIFin) / (1. - np.exp(-KNitro * LAIFin)) # SLNB calculated from exponential N profile in canopy
    SLNNTop = (NitLVI + 0.001 * NitLVI) * KNitro / (1. - np.exp(-KNitro * LAIFin)) # Value of SLNT with small plant-N increment
    RSLNBtm = (SLNBtmC - SLNBIni) / GrowP  
    
    # FVPD
    FVPDAdj = Functions.INSWTICH(C3C4Ty, 0.195127, 0.116214) # Slope for linear effect of VPDL on Ci/Ca (VPDL: Air-to-leaf vapour pressure deficit)
  
    # Canopy Photosynthesis, Transpiration, and Soil Evaporation 
    [PPHOTOCAN, APHOTOCANS, APHOTOCAN_N, APHOTOCAN, PTranspCan, ATranspCan, PEvapSoil, AEvapSoil, DiffS, DiffSU, DiffSH, PARD] = Functions.CANOPY_PHOTO_TRANSPIRATION(StressF, SunCont, SinLatD, CosLatD, DayLen, DSinB, DayDegT, TmaxD, TminD, VapPres, WindSp, C3C4Ty, LAIIni, TotalLAI, ElevM, LeafWd, RootD, SeedD, SoilSS, BulkD, KNitro, KWater, SLNitro, SLNTop, SLNNTop, MinSLN, DailyWSupply, CO2Atm, LodgeS, MaxAJ, ValN, ValJ, WatUse, WCUpperL, FVPDAdj)

    # Developmental Stage & Cumulative Thermal Units  
    TDU = Functions.THERMAL_TIME(DayShift, TmaxD, TminD, max(0., DiffS), DayLen, BaseTD, OptTD, CritTD, SenTD)
    DevRate = Functions.PHENOLOGY(DayShift, SlopeV, DayDegL, StartP, EndP, SenSc, DevDurV, DevDurR, TDU);



    
    
    # Biomass calculations
    CropWS = CropSC / CFormO
    LeafW = LeafCIni / CarFV
    StemW = StemSC / CarFV + RootVC / 0.444
    RootW = RootSC / CarFV + RootVR / 0.444
    ShootW = LeafW + StemW + CropWS
    TotalW = ShootW + RootW
    LeafWD = LeafCD / CarFV
    ShootWH = ShootW + (LeafWD - LeafCDS / CarFV)
    HarvI = CropWS / ShootWH
    
    RootWD = RootCD / CarFV
    
    # Nitrogen accumulation
    ShootN = StemN + NitLVI + CropNS
    ShootNH = ShootN + (LeafWD - LeafCDS / CarFV) * LeafNC
    TotalN = ShootN + NitRT
    
    # Maintenance and respiration (g CO2 m-2 d-1)
    RespM = max(min(44./12.*0.218*(TotalN-ShootW*LeafNC-RootW*MinRCN), APHOTOCAN-1.E-5-RspMul), 0.) #residual maintenance,
    RespT = max(0., min(APHOTOCAN-1.E-5,RspMul) + RespM) #Non-growth components of respiration, excluding the cost of N fixation
    
    # Placeholder for nitrogen fixation
    NitFix = 0
    
    RespX = 44. / 12. * (CarbonF * NitFix)
    
    # Current Assimilates for growth
    CurAssim = APHOTOCAN - RespT - RespX
    
    # Nitrogen concentration in biomass
    LeafNCB = NitLVI / LeafW  # Nitrogen concentration in leaf
    ShootNC = ShootN / ShootW
    RootNC = NitRT / RootW
    OtherNC = Functions.INSWTICH(-CropWS, CropNS / max(Functions.AVOID_ZERO_DIVISION(CropWS), 0), 0)
    PlantNC = TotalN / TotalW
    
    # Estimation of seed properties
    EndSD = Functions.INSWTICH(GrowP, DrySI, 1.) # Development stage for end of seed-number determining period
    NitRes = NitREndF + (NitREndE - NitREndF) * (EndSD - 1.0) / Functions.AVOID_ZERO_DIVISION(min(DayShift, EndSD) - 1) 
    TotalSN = NitRes / PrePN / SeedNC / SeedWt * 1000 # Thousand-seed weight
    
    # Leaf senescence
    LeafLVM = (LAIIni - min(LAIIni, LAIAdj)) / LeafSA / GrowP
    LeafLV = min (LeafW-1.E-5, LeafLVM+Functions.INSWTICH(EndSD-DayShift, LeafLVM) *0.03*LeafW)
    LeafNLev = min(LeafLV, LeafLVM) * LeafNC + (LeafLV-min(LeafLV, LeafLVM)) * LeafNCB
    LeafCLV = LeafLV*CarFV
    
    # Amount of seed protein
    SeedP = 6.25 * CropWS * OtherNC
    
    # Carbon accumulation
    CarbShootH = LeafCIni + StemSC + RootVC + CropSC
    CarbRootT = RootSC + RootVR
    CarbTotal = CarbShootH + CarbRootT
    
    # Crop nitrogen demand
    ShootSA = 12.0 / 44.0 * YielGV * (APHOTOCAN - RespT - RespX) / CarbShootH
    RespMN = max(0., min(APHOTOCAN-1.E-5,RspMul) + max(min(44./12.*0.218*(1.001*TotalN-ShootW*LeafNC-RootW*MinRCN), APHOTOCAN-1.E-5-RspMul),0.))
    ShootSAN = 12.0 / 44.0 * YielGV * (APHOTOCAN - RespMN - RespX) / CarbShootH
    NitDemR = max(0.0, (ShootSAN - ShootSA) / (0.001 * TotalN / CarbTotal))
    NitDemA = CarbRootT * ShootSA**2/Functions.AVOID_ZERO_DIVISION(NitDemR)
    
    LeafNCR = LeafNC * np.exp(-0.4 * DayShift)
    NitDemD = Functions.INSWTICH(DayShift - 1.0, ShootW * (LeafNCR - ShootNC) * (1.0 + NitRT / ShootN) / GrowP, 0.0)
    NitDemAdj = Functions.INSWTICH(LeafNCB - 1.5 * LeafNC, max(NitDemA, NitDemD), 0)
    NitDemand = Functions.INSWTICH(MinSLN - SLNitro + 1.0e-5, min(NUpTk, NitDemAdj), 0) # Crop nitrogen demand
    
    # Nitrogen partitioning
    NitCR = Functions.INSWTICH(SLNTop - MinSLN, 0, min(NUpTk, NitDemA)) / ( YielGV * (APHOTOCANS - RespT - RespX) * 12 / 44 )
    FNShoot = 1 / (1 + NitCR * NitDemR / ShootSA * CarbShootH / CarbRootT * NitRT / ShootN)
    
    # Carbon supply for shoot & root growth
    FCSHoot = 1 / (1 + NitCR * NitDemR / ShootSA) # Fraction of new carbon partitioned to shoot
    DCShootS = 12.0 / 44.0 * FCSHoot * CurAssim # Daily carbon supply from current photosynthesis for shoot growth
    DCRootS = 12 / 44 * (1 - FCSHoot) * CurAssim # Daily carbon supply from current photosynthesis for root growth
    
    # Daily carbon flow for seed filling
    FDS = Functions.SINK_GROWTH(DevRate, 1., MaxES * 1., Functions.LIMIT_FUNCTION(1., 2., DayShift) - 1)
    [SeedCDec, SeedCSupply, FlowCropS] = Functions.SINK_SOURCE(DayShift, 1., TotalSN * SeedWt * CFormO, YieldGO, FDS, DCDSR, DCShootS, GrowP)

    
    # Daily carbon flow for stem growth
    DCStemT = DCShootS - FlowCropS # Daily carbon supply from current photosynthesis for structural stem growth
    IFStemH = Functions.LIMIT_FUNCTION(0, 1, DCStemT / Functions.AVOID_ZERO_DIVISION(DCDTP)) # Integral factor of stresses on plant height growth
    FDHt = Functions.SINK_GROWTH(DevRate, (1. + DrySI) / 2., MaxEH * (1. + DrySI) / 2., min((1. + DrySI) / 2., DayShift))
    [SeedCTec, SeedCSupplyT, FlowCropST] = Functions.SINK_SOURCE(DayShift, 0, CropHt * MaxHt * CarFV, YielGV, FDHt * IFStemH, DCDTR, DCStemT, GrowP)

    
    PlantHtR = min(MaxHt - ElevM, FDHt * MaxHt * IFStemH) # rate of HT (ElevM: Plant height)
    RCarbonDTP = (SeedCTec - DCDTP) / GrowP  # Carbon demand for structural stem growth at the previous time step
    
    # Root senescence
    RootKCN = -np.log(0.05) / 6.3424 / CarFV / RootB / RootDm ## Extinction coefficient of root nitrogen 
    RootSCN = 1 / RootKCN * np.log(1.0 + RootKCN * max(0.0, (NitRT * CarFV - RootVR * MinRCN)) / MinRCN)
    RootLCT = max(min(RootSC - 1.0e-4, RootSC - min(RootSCN, RootSC)), 0.0) / GrowP
    RootLWT = RootLCT / CarFV
    RootLNT = RootLWT * MinRCN
    
    # Carbon partitioning
    FCropS = FlowCropS / DCShootS
    FStemT = Functions.INSWTICH(DayShift - (EndSD + 0.2), FlowCropST / DCShootS, 0)
    FLeafV = (1.0 if (LAIAdj - LAIIni > 0 and EndSD - DayShift > 0) else 0.0) * (1.0 - FCropS - FStemT)
    
    FRootVC = 1.0 - FLeafV - FCropS - FStemT
    FRootVR = Functions.INSWTICH(RootSCN - RootSC, 1.0, 0.0)
    
    # Dynamics of carbon-reserve pool
    GapCR = max(0, SeedCSupply - DCShootS)
    
    CREMSI = min(0.94 * RootVC, RootVC / Functions.AVOID_ZERO_DIVISION(RootVC + RootVR) * GapCR) / 0.94
    CREMRI = min(0.94 * RootVR, RootVR / Functions.AVOID_ZERO_DIVISION(RootVC + RootVR) * GapCR) / 0.94
    CREMS = Functions.INSWTICH(SeedCSupply - DCShootS, 0, CREMSI)
    CREMR = Functions.INSWTICH(SeedCSupply - DCShootS, 0, CREMRI)
    
    RRootVC = FRootVC * DCShootS - CREMS
    RRootVR = FRootVR * DCRootS - CREMR
    
    # Carbon production rate
    RCRootST = 12 / 44 * CurAssim * (1 - FCSHoot) * (1 - FRootVR) * YielGV - RootLCT
    RCStemST = 12.0 / 44.0 * CurAssim * FCSHoot * FStemT * YielGV
    RCCropS = 12.0 / 44.0 * CurAssim * FCSHoot * FCropS * YieldGO + 0.94 * (CREMS + CREMR) * YieldGO # Rate of change in Crop Storage
    RCLeafV = 12.0 / 44.0 * CurAssim * FCSHoot * FLeafV * YielGV - LeafCLV
    
    # Biomass (2)
    RWeightRoot = RCRootST / CarFV + RRootVR / 0.444
    RWeightCropS = RCCropS / CFormO
    RWeightLeaf = RCLeafV / CarFV
    RWeightStem = RCStemST / CarFV + RRootVC / 0.444
    
    # (RATES) Daily carbon flow for stem growth and seed filling
    RCDSR = max(0.0, (SeedCDec - RCCropS / YieldGO)) - (FlowCropS - min(SeedCDec, DCShootS))
    RCDTR = max(0., (SeedCTec - RCStemST / YielGV)) - (FlowCropST - min(SeedCTec, DCStemT))
    
    # Soil organic carbon dynamics
    TempFact = 47.9 / (1. + np.exp(106. / (SoilTIni + 18.3)))
    MoistFact = Functions.LIMIT_FUNCTION(0.2, 1.0, 0.2 + 0.8 * (WUpLI + WLoLI) / 10. / 150. / (FieldC - MinWC))
    DPMRConst = Functions.INSWTICH(NNUpLIni + NUpLIni + NNLoLIni + NLoLIni - AirR - NRatio, 0., DecPM)
    RPMRConst = Functions.INSWTICH(NNUpLIni + NUpLIni + NNLoLIni + NLoLIni - AirR - NRatio, 0., RootPM)
    
    CBioH = 1.67 * (1.85 + 1.60 * np.exp(-0.0786 * SoilTy))
    
    CNDPM = (DecPI + RootPMI) / Functions.AVOID_ZERO_DIVISION(DPNitI + RPNitI)
    
    DecBio = BioIni * (1 - np.exp(-TempFact * MoistFact * BioOM / 365)) / TillCP
    DecHum = HumIni * (1 - np.exp(-TempFact * MoistFact * HumusR / 365)) / TillCP
    
    DPMRate = Functions.INSWTICH(1.0 / Functions.AVOID_ZERO_DIVISION(CNDPM) - 1.0 / (8.5 * (1.0 + CBioH)), DPMRConst, DecPM)
    RPMRate = Functions.INSWTICH(1.0 / Functions.AVOID_ZERO_DIVISION(CNDPM) - 1.0 / (8.5 * (1.0 + CBioH)), RPMRConst, RootPM)
    DecDPM = DecPI * (1.0 - np.exp(-TempFact * MoistFact * DPMRate / 365)) / TillCP
    DecRPM = RootPMI * (1.0 - np.exp(-TempFact * MoistFact * RPMRate / 365)) / TillCP
    DecDPNit = DPNitI * (1.0 - np.exp(-TempFact * MoistFact * DPMRate / 365)) / TillCP
    DecRPNit = RPNitI * (1.0 - np.exp(-TempFact * MoistFact * RPMRate / 365)) / TillCP
    
    MDNit = 1.0 / 8.5 * (DecBio + DecHum) + DecDPNit + DecRPNit - 1.0 / 8.5 / (1.0 + CBioH) * (DecDPM + DecRPM + DecBio + DecHum)
    MDNUpL = (1. - np.exp(-0.065 * RootD)) * MDNit
    MDNLoL = np.exp(-0.065 * RootD) * MDNit
    
    # Soil organic nitrogen dynamics
    MinUpL = Functions.INSWTICH(MDNit, -min((NUpLIni - RootD / 150. * AirR) / TillCP, -MDNUpL), MDNUpL)
    MinLoL = Functions.INSWTICH(MDNit, -min((NLoLIni - (150. - RootD) / 150. * AirR) / TillCP, -MDNLoL), MDNLoL)
    
    MinNUpL = Functions.INSWTICH(MDNit, -min(NNUpLIni / TillCP, -MDNUpL + MinUpL), 0.)
    MinNLoL = Functions.INSWTICH(MDNit, -min(NNLoLIni / TillCP, -MDNLoL + MinLoL), 0.)
    
    MoistUpL = Functions.LIMIT_FUNCTION(0.2, 1.0, 0.2 + 0.8 * WUpLI / 10. / RootD / (FieldC - MinWC))
    MoistLoL = Functions.LIMIT_FUNCTION(0.2, 1.0, 0.2 + 0.8 * WLoLI / 10. / (150. - RootD) / (FieldC - MinWC))
    
    NitrUpL = max(0., (NUpLIni + MinUpL * TillCP - RootD / 150 * AirR)) * (1 - np.exp(-TempFact * MoistUpL * 0.6 / 7)) / TillCP
    NitrLoL = max(0., (NLoLIni + MinLoL * TillCP - (150 - RootD) / 150 * AirR)) * (1 - np.exp(-TempFact * MoistLoL * 0.6 / 7)) / TillCP
    
    ResCO2 = CBioH / (1.0 + CBioH) * (DecDPM + DecRPM + DecBio + DecHum)
    
    DeniUpL = .0005 * max(0., NNUpLIni + MinNUpL * TillCP - RootD / 150. * NRatio) * ResCO2 * (1. - np.exp(-0.065 * RootD))
    DeniLoL = .0005 * max(0., NNLoLIni + MinNLoL * TillCP - (150. - RootD) / 150. * NRatio) * ResCO2 * np.exp(-0.065 * RootD)
    
    FWStress = min(1, WUpLI / (RootD * 10 * (FieldC - MinWC)))
    
    # Total nitrogen in the soil
    TotalNA = NUpLIni + NLoLIni
    TotalNN = NNUpLIni + NNLoLIni
    NMiner = TotalNA + TotalNN
    Volatil = Functions.INSWTICH(Precip - 1., 0.15, 0.) * StressF
    
    # Crop nitrogen uptake
    NSupplyA = max(0., NUpLIni + (MinUpL - NitrUpL) * TillCP - RootD / 150. * AirR) / TillCP
    NSupplyN = max(0., NNUpLIni + (MinNUpL - DeniUpL) * TillCP - RootD / 150. * NRatio) / TillCP * FWStress
    NSupplyApp = Functions.INSWTICH(WetIdx, NitApp, NSupplyA)
    NSupplyNit = Functions.INSWTICH(WetIdx, NitFix, NSupplyN)
    NSupplyTotal = NSupplyApp + NSupplyNit
    
    NUptakeA = min(NSupplyA, NSupplyA / max(1e-10, NSupplyTotal) * max(0, NitDemand / TillCP))
    NUptakeN = min(NSupplyN, NSupplyN / max(1e-10, NSupplyTotal) * max(0, NitDemand / TillCP))
    NUptakeTotal = max(0, NUptakeA + NUptakeN + min(NitDemand / TillCP))
    
    RDemandP = (NitDemand - NitDemP) / GrowP
    RNSupplyP = (NSupplyTotal - NSupplyPP) / GrowP
    
    # Estimation of seed properties (Rates)
    RNRes = NUptakeTotal - (LeafNC * (RCLeafV + LeafCLV) + MinRCN * (RCRootST + RootLCT) + StemNC * RCStemST) / CarFV
    RNResEndE = Functions.INSWTICH(DayShift - EndSD, RNRes, 0.0)
    RNResEndF = Functions.INSWTICH(DayShift - 1.0, RNRes, 0.0)
    
    # Nitrogen accumulation (Rates)
    [RCNitRT, RCNitST, RCNitLV, RTCNitLV, RCNitSO] = Functions.NITROGEN_DYNAMIC(FNShoot, NUptakeTotal, RWeightStem, StemNC, LeafNC, MinRCN, LeafNCB, RootNC, NitLVI, NitRT, LeafW, RootW, GrowP, CarbonB, CarbonX, TempM, DayShift, SeedNC, RWeightCropS, LeafNLev, RootLNT)
    
    # Leaf area development (Rates)
    RLAI = Functions.LAI_RATE(DayShift, LeafSA, RWeightLeaf, LAIIni, KNitro, NitLVI, RCNitLV, SLNBIni, RSLNBtm)
    
    # Rooting depth (Rates)
    RootKRate = -np.log(0.05) / RootDm # Extinction coefficient of root weight density over the soil depth
    RootDepthRate = Functions.INSWTICH(RootD - RootDm, min((RootDm - RootD) / GrowP, (RWeightRoot + RootLWT) / (RootB + RootKRate * (RootW + RootWD))), 0) # Rooting depth rate
    
    
        # Maintenance and respiration adjustments
    ResMD = 0.06 * (1 - FCSHoot) * CurAssim
    ResMS = 0.06 * 0.05 / 0.454 * YielGV * CurAssim
    ResMN = 44 / 12 * 2.05 * NUptakeN
    ResMA = 44.0 / 12.0 * 0.17 * NUptakeA
    RateRspMul = (ResMN + ResMA + ResMS + ResMD - RspMul) / GrowP
    
    GrowthR = 44.0 / 12.0 * (((1.0 - YielGV) / YielGV * (RCLeafV + RCStemST + RCRootST + LeafCLV + RootLCT)) + ((1.0 - YieldGO) / YieldGO * RCCropS)) # Daily growth respiration
    
    TotalResp = RespT + RespX + GrowthR + 44.0 / 12.0 * 0.06 * (CREMS + CREMR) # Daily total respiratory cost 
    
    # Soil temperature dynamics
    AvgSoilT = ((AvgTemp + DiffS) + NAvgTemp) / 2.
    RateSoilT = (AvgSoilT - SoilTIni) / TillCT
    
    # Daily and total C and N returns from crop to soil
    LeafDecS = (LeafCD - LeafCDS) / 10. * (AvgSoilT - BaseTD) / (OptTD - BaseTD)
    LitCarbon = RootLCT + LeafDecS
    LitNitro = RootLNT + LeafDecS / CarFV * LeafNC * PanLS
    NitroRetS = LitNitroT + Functions.INSWTICH(DayShift - 2, 0, NitLVI + StemN + NitRT + (LeafCD - LeafCDS) / CarFV * LeafNC * (1 + PanLS) / 2)
    
    # Carbon balance check for the crop
    CarbonCheckIn = CarbTotal + LeafCD + RootCD - LeafCIni - RootCIni
    CarbonCheck = (CarbonCheckIn - (CurAssim - TotalResp) * 12.0 / 44.0) / Functions.AVOID_ZERO_DIVISION(CarbonCheckIn) * 100
    
    # Nitrogen balance check for the crop
    NitroCheckIn = TotalN + LeafND + RootND - NitLVI - NitRTI
    NitroCheck = (NitroCheckIn - NUptakeTotal) / Functions.AVOID_ZERO_DIVISION(NUptakeTotal) * 100
    
    # Soil water availability and water balance adjustments
    Irrig = 0
    RainFir = Precip + Irrig
    RootZoneUpL = min(10.0 * (MaxWC - WCUpperL) * RootD / TillCP, RainFir)
    RootZoneLoL = min(10. * (MaxWC - WCLowerL) * (150. - RootD) / TillCP, RainFir - RootZoneUpL)
    WaterUptakeG = max(0., RainFir - RootZoneUpL - RootZoneLoL)
    WaterLoL = RootZoneLoL - 10.0 * (WCLowerL - MinWC) * RootDepthRate
    WaterUpL = RootZoneUpL + 10.0 * (WCLowerL - MinWC) * RootDepthRate - Functions.INSWTICH(WetIdx, 0.0, ATranspCan + AEvapSoil) + 0.1
    
    # Soil organic nitrogen adjustments
    FernA = 0
    FernN = 0
    
    RateSoilFernA = FernA - SoilFernA / 3
    LeachAUL = max(0., (NSupplyNit - NUptakeN) * TillCP - RootD /150.*NRatio) * min((RainFir - RootZoneUpL) / MaxWC / RootD / 10., 1.)
    LayerNitA = RootDepthRate / (150.0 - RootD) * NLoLIni
    LayerNitN = RootDepthRate / (150.0 - RootD) * NNLoLIni
    
    # Soil organic carbon rates
    RateDPM = LitCarbon * DecPM / (1. + DecPM) - DecDPM
    RateRPM = LitCarbon * 1. / (1. + DecPM) - DecRPM
    RateDPNit = LitNitro / (1. + 40. / DecPM / 100.) - DecDPNit
    RateRPNit = LitNitro / (1. + 100. * DecPM / 40.) - DecRPNit
    RateBio = 0.46 / (1.0 + CBioH) * (DecDPM + DecRPM + DecBio + DecHum) - DecBio
    RateHum = 0.54 / (1.0 + CBioH) * (DecDPM + DecRPM + DecBio + DecHum) - DecHum
    RateNUpL = FernA + MinUpL + LayerNitA - Functions.INSWTICH(WetIdx, 0.0, NUptakeA) - NitrUpL - Volatil
    RateNLoL = MinLoL - LayerNitA - NitrLoL
    RateNNLoL = LeachAUL + MinNLoL + NitrLoL - LayerNitN - DeniLoL 
    RateNNUpL = FernN + MinNUpL + NitrUpL + LayerNitN - Functions.INSWTICH(WetIdx, 0.0, NUptakeN) - DeniUpL - LeachAUL
    
    # Incrementing state variables with their respective rates
    DayShift += DevRate
    CumTDU += TDU
    LeafCIni += RCLeafV
    LeafCD += LeafCLV
    StemSC += RCStemST
    CropSC += RCCropS
    RootSC = RCRootST + RootCIni
    RootCD += RootLCT
    LeafCDS += LeafDecS
    NitRT += RCNitRT
    StemN += RCNitST
    NitLVI += RCNitLV
    CropNS += RCNitSO
    TotalNLV += RTCNitLV
    LeafND += LeafNLev
    RootND += RootLNT
    RootVC += RRootVC
    RootVR += RRootVR
    NitREndE += RNResEndE
    NitREndF += RNResEndF
    DCDSR += RCDSR
    DCDTR += RCDTR
    SLNBIni += RSLNBtm
    LAIIni += RLAI
    RspMul += RateRspMul
    NitDemP += RDemandP
    NSupplyPP += RNSupplyP
    DCDTP += RCarbonDTP
    ElevM += PlantHtR
    RootD += RootDepthRate
    TotalDAPAR += PARD
    TotalAPCan += APHOTOCAN
    TotalResp += TotalResp
    TotalTranspCan += ATranspCan
    TotalNUptake += NUptakeTotal
    LitNitroT += LitNitro
    SoilTIni += RateSoilT
    WUpLI += WaterUpL
    WLoLI += WaterLoL
    DecPI += RateDPM
    RootPMI += RateRPM
    BioIni += RateBio
    HumIni += RateHum
    DPNitI += RateDPNit
    RPNitI += RateRPNit
    TotalNNLoL += LeachAUL
    NUpLIni += RateNUpL
    NLoLIni += RateNLoL
    NNUpLIni += RateNNUpL
    NNLoLIni += RateNNLoL
    SoilFernA += RateSoilFernA

    
    
    # TABLE_OUT[I,0]=TIME
    # TABLE_OUT[I,1]=DCDSC
    # TABLE_OUT[I,2]=DCDS
    # TABLE_OUT[I,3]=FLWCS
    # TABLE_OUT[I,4]=DCDT
    # TABLE_OUT[I,5]=FLWCT
    # TABLE_OUT[I,6]=RLAI
    # TABLE_OUT[I,7]=LAI
    


