import math
from equations import meteo
from equations import globals

# header file: solar.h


# Extinction coefficient for direct fluxes. th is solar time
def Kdr(th):
    dElev = meteo.SolarElevation(th)
    if dElev < 0.00000001:
        return 0.0 # sun is below horizon
    return (0.5 / math.sin(dElev))

# Extinction coefficient for diffuse fluxes
def Kdf():
    dL = math.sqrt(globals.PARAM.lai)
    return ((1.0 + 0.1174 * dL) / (1.0 + 0.3732 * dL))

# Hourly intercepted of total solar radiation (W m-2)
#   th (solar time),
#   dfrad (diffuse irradiance), drrad (direct irradiance)
#   where both dfrad and drrad will be set inside this function
def InterceptHourRad(th, dfrad, drrad):
    P = 0.11 # reflection coefficient
    S = math.sqrt(0.50) # scatter correction
    dDfRad = 0.0
    dDrRad = 0.0
    dDfRad, dDrRad=meteo.HourRad(th, dDfRad, dDrRad)
    drrad = (1.0 - P) * dDrRad * (1.0 - math.exp(-Kdr(th) * S * globals.PARAM.lai))
    dfrad = (1.0 - P) * dDfRad * (1.0 - math.exp(-Kdf() * S * globals.PARAM.lai))
    return(dfrad, drrad)
# Daily intercepted of total solar radiation (J m-2 day-1)
#   dfrad (diffuse irradiance), drrad (direct irradiance)
#   where both dfrad and drrad will be set inside this function

def InterceptDayRad(dfrad, drrad):
    print("**** HERE needs an edit .... solar.py .... *********************************************")
    # numerical integration using 5-point Gaussian method
    ABS = [0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899]
    WGT = [0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635]
    dTsr = meteo.Sunrise()
    dIval = meteo.Sunset() - dTsr
    dfra = 0.0
    drrad = 0.0
    for i in range(0, 5):
        dHr = dTsr + ABS[i] * dIval
        dDfrad = 0.0
        dDrrad = 0.0
        dDfrad, dDrrad=InterceptHourRad(dHr, dDfrad, dDrrad)
        # x 3600 to convert sec to hour:
        dfrad = dfrad + dDfrad * 3600 * WGT[i] * dIval
        drrad = drrad + dDrrad * 3600 * WGT[i] * dIval
    return (dfrad, drrad)


# Absorbed PAR (W m-2) by sunlit and shaded leaves.
#   th (solar time), 
#   both sunlit and shaded will be set inside this function
def AbsorbedHourPAR(th, sunlit, shaded):
    P = 0.04 # reflection coefficient
    A = 0.80 # scatter coefficient
    S = math.sqrt(A) # scatter correction
    dDfPar = 0.0
    dDrPar = 0.0
    dDfPar, dDrPar=meteo.HourRad(th, dDfPar, dDrPar)

    dDfPar = 0.5 * dDfPar # 50% total radiation is PAR
    dDrPar = 0.5 * dDrPar

    dKdr = Kdr(th)
    dI = (1.0 - P) * dDrPar
    dIpdr = dI * math.exp(-dKdr * S * globals.PARAM.lai) # total direct
    dIpdrdr = dI * math.exp(-dKdr * globals.PARAM.lai) # direct of direct
    dIpdra = (dIpdr - dIpdrdr) / 2 # scatter beams
    dN = Kdf() * S * globals.PARAM.lai
    dI = (1.0 - P) * dDfPar
    dIpdf = dI * (1.0 - math.exp(-dN)) / dN # diffuse

    sunlit = A * (dKdr * dDrPar + dIpdf + dIpdra)
    shaded = A * (dIpdf + dIpdra)
    return (sunlit, shaded)

# Leaf area index for sunlit and shaded leaves.
#   th (solar time).
#   Both sunlit and shaded will be set inside this function
def LAI(th, sunlit, shaded):
    sunlit = 0.0
    dKdr = Kdr(th)
    if dKdr > 0.0:
        sunlit = (1.0 - math.exp(-dKdr * globals.PARAM.lai)) / dKdr
    shaded = globals.PARAM.lai - sunlit
    return (sunlit, shaded)


