import math
from equations import solar
from equations import meteo
from equations import et

# header file: assim.h



# Changes the given model parameter based on current temperature.
#   k25 (model parameter at 25 deg. C), q10 (Q10 value),
#   lftemp (leaf temperature, deg. C)
def Q10(k25, q10, lftemp):
    return (k25 * q10 ** ((lftemp - 25.0) / 10.0))

# Rubisco capacity rate (umol m-2 s-1)
#   lftemp (leaf temperature, deg. C)
def Vcmax(lftemp):
    VCMAX = 200.0 # max rate
    dExp = 1.0 + math.exp(0.128 * (lftemp - 40.0))
    return (Q10(VCMAX, 2.4, lftemp) / dExp)

# CO2 compensation point (umol mol-1)
#   lftemp (leaf temperature, deg. C)
def CO2CompensationPt(lftemp):
    TAU = 2600.0 # specificity factor
    OA = 210000.0 # O2 concentration
    dTau = Q10(TAU, 0.57, lftemp)
    return (OA * 0.5 / dTau)

# Gross photosynthesis limited by light (umol m-2 s-1)
#   par (PAR, umol m-2 s-1), lftemp (leaf temperature, deg. C),
#   ci (intercellular CO2 concentration, umol mol-1)
def LightLimited(par, lftemp, ci):
    DQE = 0.06 # quantum yield
    ALPHA = 0.8 # absorption fraction
    dT = CO2CompensationPt(lftemp)
    dA = ALPHA * DQE * par * (ci - dT)
    dB = ci + 2.0 * dT
    return (dA / dB)

# Gross photosynthesis limited by Rubisco capacity (umol m-2 s-1)
#   lftemp (leaf temperature, deg. C),
#   ci (intercellular CO2 concentration, umol mol-1)
def RubiscoLimited(lftemp, ci):
    KCMAX = 300.0 # O2
    KOMAX = 300000.0 # CO2
    OA = 210000.0 # O2 concentration
    dKc = Q10(KCMAX, 2.1, lftemp)
    dKo = Q10(KOMAX, 1.2, lftemp)
    dKm = dKc * (1.0 + OA / dKo)
    dA = Vcmax(lftemp) * (ci - CO2CompensationPt(lftemp))
    dB = ci + dKm
    return (dA / dB)

# Gross photosynthesis limited by sucrose sink (umol m-2 s-1)
#   lftemp (leaf temperature, deg. C)
def SinkLimited(lftemp):
    return (0.5 * Vcmax(lftemp))

# Gross leaf photosynthesis (umol m-2 s-1)
#   par (PAR, umol m-2 s-1), lftemp (leaf temperature, deg. C),
def LeafAssim(par, lftemp):
    CI = 245.0 # internal CO2 concentration (umol mol-1)
    dMin = RubiscoLimited(lftemp, CI)
    dJc = LightLimited(par, lftemp, CI)
    dJs = SinkLimited(lftemp)

    # find the most limiting factor to assimilation
    if dMin > dJc:
        dMin = dJc

    if dMin > dJs:
        dMin = dJs

    return dMin

# Gross canopy photosynthesis (per unit ground area)
#   (umol CO2 m-2 s-1).
#   th (local solar time, hour), lai (leaf area index, m2 m-2),
#   lftemp (leaf temperature, deg. C)
def HourCanopyAssim(th, lftemp):
    # absorbed PAR
    dQsl = 0.0
    dQsh = 0.0
    dQsl, dQsh=solar.AbsorbedHourPAR(th, dQsl, dQsh)
    dQsl = 4.55 * dQsl # convert to umol m-2 s-1
    dQsh = 4.55 * dQsh
    # leaf assimilation:
    dAsl = LeafAssim(dQsl, lftemp)
    dAsh = LeafAssim(dQsh, lftemp)
    # sunlit and shaded LAI
    dLsl = 0.0
    dLsh = 0.0
    dLsl, dLsh=solar.LAI(th, dLsl, dLsh)
    return (dAsl * dLsl + dAsh * dLsh) # in umol CO2 m-2 s-1

# Daily gross canopy photosynthesis (per unit ground area)
#   (umol CO2 m-2 day-1)
def DayCanopyAssim():
    # 5-point gauss integration over sunrise to sunset:
    ABS = [0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899]
    WGT = [0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635]

    dTsr = meteo.Sunrise()
    dIval = meteo.Sunset() - dTsr
    dAssim = 0.0
    for i in range(0, 5):
        dHour = dTsr + ABS[i] * dIval # current hour
        dET = 0.0
        dETs = 0.0
        dETc = 0.0
        dH = 0.0
        dHs = 0.0
        dHc = 0.0
        dET, dETs, dETc, dH, dHs, dHc=et.HourHeatFluxes(dHour, dET, dETs, dETc, dH, dHs, dHc)
        dTf = et.CanopyTemp(dHour, dH, dHc)
        dAn = HourCanopyAssim(dHour, dTf)
        # x 3600 to convert sec to hour:
        dAssim = dAssim + dAn * 3600.0 * WGT[i] * dIval

    return dAssim


