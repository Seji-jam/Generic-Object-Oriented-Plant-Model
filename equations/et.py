import math
from equations import meteo
from equations import globals
from equations import solar

# header file: et.h


# Slope of saturated vapor pressure with temperature (mbar/K).
#   temp (temperature)
def SlopeSvp(temp):
    dA = 25029.4 * math.exp(17.269 * temp / (temp + 237.3))
    dB = (temp + 237.3) ** 2
    return (dA / dB)

# Vapor pressure deficit at reference height (mbar).
#   temp (temperature)
def VpdRef(temp):
    return (meteo.Svp(temp) - meteo.HourVP(temp))

# Vapor pressure deficit at mean canopy flow (mbar).
#   a (total energy available), et (total ET),
#   raa (aerodynamic resistance between mean canopy flow
#   and reference height), temp (temperature)
def VpdMcf(a, et, raa, temp):
    dDelta = SlopeSvp(temp)
    dA = (raa / 1221.09) * (dDelta * a - (dDelta + 0.658) * et)
    return (VpdRef(temp) + dA)

# Sensible heat for soil (W m-2).
#   rsa (aerodynamic resistance between mean canopy flow
#   and soil), rss (soil resistance), vpdmcf (vapor pressure
#   deficit at mean canopy flow), as (energy available to soil),
#   slopesvp (slope of saturated vapor pressure
#   with temperature)
def Hs(rsa, rss, vpdmcf, as_, slopesvp):
    dA = 0.658 * as_ * (rss + rsa) - 1221.09 * vpdmcf
    dB = slopesvp * rsa + 0.658 * (rss + rsa)
    return (dA / dB)

# Sensible heat for crop (W m-2).
#   rca (aerodynamic resistance between mean canopy flow
#   and canopy), rcs (canopy resistance), vpdmcf (vapor pressure
#   deficit at mean canopy flow), ac (energy available to
#   canopy), slopesvp (slope of saturated vapor pressure
#   with temperature)
def Hc(rca, rcs, vpdmcf, ac, slopesvp):
    dA = 0.658 * ac * (rcs + rca) - 1221.09 * vpdmcf
    dB = slopesvp * rca + 0.658 * (rcs + rca)
    return (dA / dB)

# Latent heat for soil (W m-2).
#   rsa (aerodynamic resistance between mean canopy flow
#   and soil), rss (soil resistance), vpdmcf (vapor pressure
#   deficit at mean canopy flow), as (energy available to soil),
#   slopesvp (slope of saturated vapor pressure
#   with temperature)
def ETs(rsa, rss, vpdmcf, as_, slopesvp):
    dA = slopesvp * as_ + (1221.09 * vpdmcf / rsa)
    dB = slopesvp + 0.658 * (rss + rsa) / rsa
    return (dA / dB)

# Latent heat for crop (W m-2).
#   rca (aerodynamic resistance between mean canopy flow
#   and canopy), rcs (canopy resistance), vpdmcf (vapor pressure
#   deficit at mean canopy flow), ac (energy available to
#   canopy), slopesvp (slope of saturated vapor pressure
#   with temperature)
def ETc(rca, rcs, vpdmcf, ac, slopesvp):
    dA = slopesvp * ac + (1221.09 * vpdmcf / rca)
    dB = slopesvp + 0.658 * (rcs + rca) / rca
    return (dA / dB)

# Soil heat flux (W m-2).
#   th (local solar time), rsn (solar radiation reaching the soil)
def G(th, rsn):
    dSolarInc = (globals.Globals.gPI/2.0) - meteo.SolarElevation(th)
    return (0.35 * math.cos(dSolarInc) * rsn)

# Zeroplane displacement (m)
def Zeroplane():
    return (0.64 * globals.PARAM.hgt)

# Surface roughness length (m)
def RoughLen():
    return (0.13 * globals.PARAM.hgt)

# Friction velocity (m2/s2). hgt (crop height),
#   windspd (wind speed at canopy top)
def FrictionVelocity(windspd):
    dD = Zeroplane()
    dZ = RoughLen()
    return (0.4 * windspd / math.log((globals.PARAM.hgt - dD) / dZ))

# Calculate the two aerodynamic resistances (s/m):
#   rsa (aerodynamic resistance between mean canopy flow
#   and soil), and raa (aerodynamic resistance between mean
#   canopy flow and reference height), 
#   windspd (wind speed at canopy top)
def AeroResistances(windspd, rsa, raa):
    dD = Zeroplane()
    dZ = RoughLen()
    dFricVel = FrictionVelocity(windspd)
    dKh = 0.4 * dFricVel * globals.PARAM.hgt
    dP = math.exp(-globals.PARAM.keddy * 0.004 / globals.PARAM.hgt) - math.exp(-globals.PARAM.keddy * (dZ + dD) / globals.PARAM.hgt)
    rsa = (globals.PARAM.hgt * math.exp(globals.PARAM.keddy) / (globals.PARAM.keddy * dKh)) * dP
    dP = math.exp(globals.PARAM.keddy * (1.0 - (dZ + dD) / globals.PARAM.hgt)) - 1.0
    dP = (globals.PARAM.hgt / (globals.PARAM.keddy * dKh)) * dP
    raa = math.log((globals.PARAM.refhgt - dD) / (globals.PARAM.hgt - dD)) / (0.4 * dFricVel) + dP
    return (rsa, raa)

# Boundary layer resistance (s/m), rca (resistance between
#   mean canopy flow and canopy).
#   windspd (wind speed at canopy top)
def BoundLayerResistance(windspd):
    dD = Zeroplane()
    dZ = RoughLen()
    dFricVel = FrictionVelocity(windspd)
    dUh = (dFricVel / 0.4) * math.log((globals.PARAM.hgt - dD) / dZ)
    dRca = 0.012 * globals.PARAM.lai * (1.0 - math.exp(-globals.PARAM.kwind/2.0))
    dRca = globals.PARAM.kwind / (dRca * math.sqrt(dUh / globals.PARAM.leafwidth))
    return dRca

# Canopy resistance (s/m), rcs.
#   totrad (total solar irradiance)
def CanopyResistance(totrad):
    if totrad <= 0.0:
        totrad = 0.01 # avoid dividing by zero

    dPar = 0.5 * totrad # PAR is 50% of solar rad.
    dRcs = (globals.PARAM.rstA1 + dPar) / (globals.PARAM.rstA2 * dPar)

    LAICR = 0.5 * 4.0 # critical LAI
    dLe = globals.PARAM.lai
    if globals.PARAM.lai > LAICR:
        dLe = LAICR

    return (dRcs / dLe)

# Soil resistance (s/m), rss
def SoilResistance():
    dRssdry = (2.0 * 0.02) / (globals.PARAM.porosity * 0.0000247)
    dExp = math.exp(-(1.0 / globals.PARAM.poredist) * (globals.PARAM.vwc01 / globals.PARAM.vwcsat))
    return (dRssdry * dExp)

# Net energy available to soil and canopy (W m-2).
#   th (local solar time), lai (leaf area index),
#   as (energy available to soil), ac (energy available
#   to canopy). Both as and ac will be set inside this function.
def EnergySupply(th, as_, ac):
    ac = as_ = 0.0
    dRn = meteo.HourNetRad(th)
    if dRn > 0.0:
        S = math.sqrt(0.50) # scatter correction
        dRsn = dRn * math.exp(-solar.Kdr(th) * S * globals.PARAM.lai) # at soil
        ac = dRn - dRsn
        as_ = dRn - G(th, dRsn) - ac
    return(as_, ac)

# Canopy temperature (deg. C).
#   th (local solar time),
#   h and hc - sensible heat for total and canopy
def CanopyTemp(th, h, hc):
    dWind = meteo.HourWind()
    dRsa = 0.0
    dRaa = 0.0
    dRsa, dRaa=AeroResistances(dWind, dRsa, dRaa)
    dT0 = ((h / 1221.09) * dRaa) + meteo.HourTemp(th)
    dRca = BoundLayerResistance(dWind)
    dTf = ((hc / 1221.09) * dRca) + dT0
    return dTf

# Latent and sensible heat fluxes (W m-2) at the given hour.
#   th (local solar time).
#   Variables set inside this function:
#   et, ets and etc - latent heat for total, soil and canopy
#   h, hs and hc - sensible heat for total, soil and canopy
def HourHeatFluxes(th, et, ets, etc, h, hs, hc):
    # hourly meteorological properties: /////////////////////////
    dWind = meteo.HourWind() # wind speed
    dTemp = meteo.HourTemp(th) # air temperature
    dDfrad = 0.0
    dDrrad = 0.0
    dDfrad,dDrrad= meteo.HourRad(th, dDfrad, dDrrad)
    dTotrad = dDfrad + dDrrad # solar irradiance

    # resistances: //////////////////////////////////////////////
    dRsa = 0.0
    dRaa = 0.0
#    print(globals.PARAM.hgt)
    dRsa, dRaa=AeroResistances(dWind, dRsa, dRaa)
    dRca = BoundLayerResistance(dWind)
    dRcs = CanopyResistance(dTotrad)
    dRss = SoilResistance()

#    print(dRca)


    # energy available to soil and plant: ///////////////////////
    dAs = 0.0
    dAc = 0.0
    dAs, dAc=EnergySupply(th, dAs, dAc)
#    print(dAs)
#    print(dAc)


    # calculate fluxes now: /////////////////////////////////////
    dDelta = SlopeSvp(dTemp)
    dRa = (dDelta + 0.658) * dRaa
    dRc = (dDelta + 0.658) * dRca + 0.658 * dRcs
    dRs = (dDelta + 0.658) * dRsa + 0.658 * dRss
    dCc = 1.0 + (dRc * dRa / (dRs * (dRc + dRa)))
    dCc = 1.0 / dCc
    dCs = 1.0 + (dRs * dRa / (dRc * (dRs + dRa)))
    dCs = 1.0 / dCs
    dA = dAs + dAc
    dD = VpdRef(dTemp)

    dP = dDelta * dA + ((1221.09 *dD - dDelta *dRca *dAs) / (dRaa + dRca))
    dQ = dDelta + 0.658 * (1.0 + dRcs / (dRaa + dRca))
    dPMc = dP / dQ

    dP = dDelta * dA + ((1221.09 *dD - dDelta *dRsa *dAc) / (dRaa + dRsa))
    dQ = dDelta + 0.658 * (1.0 + dRss / (dRaa + dRsa))
    dPMs = dP / dQ

    et = dCc * dPMc + dCs * dPMs

#    print(dRaa)
    dD0 = VpdMcf(dA, et, dRaa, dTemp)
    
    ets = ETs(dRsa, dRss, dD0, dAs, dDelta)

    etc = et - ets

    hs = Hs(dRsa, dRss, dD0, dAs, dDelta)
    hc = Hc(dRca, dRcs, dD0, dAc, dDelta)
    h = hs + hc
    return(et, ets, etc, h, hs, hc)

# Latent and sensible heat fluxes (J m-2 day-1) for the day.
#   Variables set inside this function:
#   et, ets and etc - latent heat for total, soil and canopy
#   h, hs and hc - sensible heat for total, soil and canopy
def DayHeatFluxes(et, ets, etc, h, hs, hc):
    # 5-point gauss integration:
    ABS = [0.0469101, 0.2307534, 0.5000000, 0.7692465, 0.9530899]
    WGT = [0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635]

    dIval = 24.0 # integration over 24 hours
    et = ets = etc = h = hs = hc = 0.0 # reset to zeros
    for i in range(0, 5):
        dHour = 0.0 + ABS[i] * dIval # current hour
        det = 0.0
        dets = 0.0
        detc = 0.0
        dh = 0.0
        dhs = 0.0
        dhc = 0.0
        det,dets, detc, dh, dhs, dhc=HourHeatFluxes(dHour, det, dets, detc, dh, dhs, dhc)


        # x 3600 to convert sec to hour:
        et = et + det * 3600.0 * WGT[i] * dIval
        ets = ets + dets * 3600.0 * WGT[i] * dIval
        etc = etc + detc * 3600.0 * WGT[i] * dIval
        h = h + dh * 3600.0 * WGT[i] * dIval
        hs = hs + dhs * 3600.0 * WGT[i] * dIval
        hc = hc + dhc * 3600.0 * WGT[i] * dIval
#    print(ets)
    return (et, ets, etc, h, hs, hc)

