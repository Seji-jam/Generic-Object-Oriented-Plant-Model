import math
from equations import globals
from equations import wthrfile

# Water content before redistribution (mm)
#   length = soil layer thickness (m)
#   wc = current water content (m3 m-3)
#   wcperc = water percolating from the above soil layer (mm day-1)
def WcBeforeRedist(length, wc, wcperc):
    dC = 1000.0 * length # for conversion
    dWcsat = globals.PARAM.vwcsat * dC
    dPercEx = 0.0 # percolation to below
    dStore = wc + wcperc # total water amount
    if dStore > dWcsat:
        dPercEx = dStore - dWcsat # total exceeds saturation
    return (wc + wcperc - dPercEx)

# Water content after redistribution (mm)
#   length = soil layer thickness (m)
#   wc = current water content (m3 m-3)
def WcAfterRedist(length, wc):
    dC = 1000.0 * length # for conversion
    dVwc = wc / dC # want unit in m3 m-3
    dExp = math.exp((globals.PARAM.alpha / globals.PARAM.vwcsat) * (globals.PARAM.vwcsat - dVwc))
    dVal = (globals.PARAM.alpha * globals.PARAM.ksat * globals.PARAM.delta) / (length * globals.PARAM.vwcsat)
    dLog = math.log(dVal + dExp)
    dWc = globals.PARAM.vwcsat - (globals.PARAM.vwcsat/globals.PARAM.alpha) * dLog
    if dWc < 0:
        dWc = 0.0 # no negative water content
    return (dWc * dC) # convert back to mm

# Actual soil evaporation (mm day-1)
#   potE = potential soil evaporation (mm day-1)
#   e01 and e02 = actual evaporation from layer 1 and 2 (mm day-1)
#                 (set within this function)
def ActualE(potE, e01, e02):
    dPow = (3.6073 * globals.PARAM.vwc01 / globals.PARAM.vwcsat) ** -9.3172
    dReduction = 1.0 / (1.0 + dPow)
    e01 = potE * dReduction * 0.26 # 26% from top layer
    e02 = potE * dReduction * 0.74 # rest from root zone
    return (e01,e02)

# Actual plant transpiration (mm day-1)
#   p = 0.5 or 0.3 for C3 and C4 plants, respectively
#   potT = potential soil transpiration (mm day-1)
#   t01 and t02 = actual transpiration from layer 1 and 2 (mm day-1)
#                 (set within this function)
def ActualT(p, potT, t01, t02):
    # critical water content
    dVwccr = globals.PARAM.vwcwp + p * (globals.PARAM.vwcsat - globals.PARAM.vwcwp)
    dReduction = (globals.PARAM.vwc02 - globals.PARAM.vwcwp) / (dVwccr - globals.PARAM.vwcwp)
    if dReduction > 1.0:
        dReduction = 1.0

    t01 = 0.0 # no active roots in top layer
    t02 = potT * dReduction # all roots in root zone
    return(t01, t02)

# Final volumetric water content (m3 m-3)
#   length = soil layer thickness (m)
#   vwc = current volumetric water content (m3 m-3)
#   ea and ta = actual evaporation and transpiration (mm day-1)
#   wcperc0 and wcperc1 = percolation from above and to below
#                         (mm day-1)
#   Note: wcperc1 is set within this function.
def Vwc(length, vwc, ea, ta, wcperc0, wcperc1):
    dC = 1000.0 * length # for conversions
    dWc = vwc * dC # want unit in mm
    dWc = WcBeforeRedist(length, dWc, wcperc0)

    dBal = WcAfterRedist(length, dWc)
    dBal = dBal - ea - ta # final water amount
    if dBal < 0.0:
        dBal = 0.0

    wcperc1 = dWc - ea - ta - dBal # percolation to below
    if wcperc1 < 0.0:
        wcperc1 = 0.0

    return (wcperc1, (dBal / dC)) # unit m3 m-3

# Final volumetric water content (m3 m-3) for the two soil
#   layers. potE and potT are the potential evaporation and
#   transpiration (mm day-1).
def DailyWaterContent(potE, potT):
    ea01 = 0.0
    ea02 = 0.0
    ea01,ea02=ActualE(potE, ea01, ea02)

    ta01 = 0.0
    ta02 = 0.0
    ta01, ta02=ActualT(0.5, potT, ta01, ta02) # assume a C3 plant

    dPerc01 = 0.0
    dPerc01,globals.PARAM.vwc01 = Vwc(globals.PARAM.len01, globals.PARAM.vwc01, ea01, ta01, wthrfile.DAILY_METEO.rain, dPerc01)
    dPerc01,globals.PARAM.vwc02 = Vwc(globals.PARAM.len02, globals.PARAM.vwc02, ea02, ta02, dPerc01, dPerc01)


# Fraction of growth reduced due to limited water. potT is potential
#   transpiration (mm day-1).
def GrowthReduction(potT):
    ta01 = 0.0
    ta02 = 0.0
    ta01, ta02=ActualT(0.5, potT, ta01, ta02) # C3 plant
    return (ta02 / potT) # fraction (0 to 1)




