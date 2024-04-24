import math
from equations import globals
from equations import wthrfile

# header file: meteo.h



# Solar declination (radians)
def SolarDeclination():
    return (-0.4093 * math.cos(2 * globals.Globals.gPI * (globals.Globals.gDoy + 10) / 365))

# Solar angle from horizontal (radians). Note: function will
#   return a -ve value if sun is below horizon. th (solar time)
def SolarElevation(th):
    dDecl = SolarDeclination()
    dA = math.sin(dDecl) * math.sin(globals.PARAM.lat)
    dB = math.cos(dDecl) * math.cos(globals.PARAM.lat)
    dHa = globals.Globals.gPI * (th - 12) / 12
    return (math.asin(dA + dB * math.cos(dHa))) # can be -ve

# Solar angle from north in an eastward direction (radians)
#   th (solar time)
def SolarAzimuth(th):
    dDecl = SolarDeclination()
    dElev = SolarElevation(th)
    dVal = (math.sin(globals.PARAM.lat) * math.sin(dElev) - math.sin(dDecl)) / (math.cos(globals.PARAM.lat) * math.cos(dElev))

    if dVal > 1.0: # problem if dVal > 1 (see *below)
        dVal = 1.0

    dAzi = math.acos(dVal) # * ensure dVal <= 1.0
    if th < 12:
        dAzi = -dAzi # -ve for before solar noon, else +ve
    return (globals.Globals.gPI + dAzi)

# Length of day (hours) - from sun rise to sun set
def Daylength():
    dDecl = SolarDeclination()
    dA = math.sin(dDecl) * math.sin(globals.PARAM.lat)
    dB = math.cos(dDecl) * math.cos(globals.PARAM.lat)
    return (24 * math.acos(-dA / dB) / globals.Globals.gPI)

# Solar time of sun set (hours)
def Sunset():
    dDecl = SolarDeclination()
    #print(dDecl)
    dA = math.sin(dDecl) * math.sin(globals.PARAM.lat)
    #print(dA)
    dB = math.cos(dDecl) * math.cos(globals.PARAM.lat)
    #print(dB)

    return (12 * math.acos(-dA / dB) / globals.Globals.gPI + 12)

# Solar time of sun rise (hours)
def Sunrise():
    return (24 - Sunset())

# Hourly extraterrestrial solar irradiance (W m-2). th (solar time)
def HourETRad(th):
    dE0 = 1 + 0.033 * math.cos(2 * globals.Globals.gPI * (globals.Globals.gDoy - 10) / 365)
    dETRad = 1370 * dE0 * math.sin(SolarElevation(th))
    if dETRad < 0.0:
        dETRad = 0.0 # no -ve radiation values
    return dETRad

# Daily extraterrestrial solar irradiance (J m-2 day-1)
def DayETRad():
    dDecl = SolarDeclination()
    dA = math.sin(dDecl) * math.sin(globals.PARAM.lat)
    dB = math.cos(dDecl) * math.cos(globals.PARAM.lat)
    dAoB = dA / dB
    dSinbe = (24.0 / globals.Globals.gPI) * (dA * math.acos(-dAoB) + dB * math.sqrt(1 - dAoB * dAoB))
    dE0 = 1.0 + 0.033 * math.cos(2 * globals.Globals.gPI * (globals.Globals.gDoy - 10) / 365)
    dScc = 1370.0 * dE0
    return (3600 * dScc * dSinbe)

# Daily solar irradiance (J m-2 day-1)
#   dfrad (diffuse irradiance), drrad (direct irradiance)
#   where both dfrad and drrad will be set inside this function
def DayRad(dfrad, drrad):
    dRelSunhr = wthrfile.DAILY_METEO.avsun / Daylength()
    dETRad = DayETRad()
    dTotrad = dETRad * (globals.PARAM.angb0 + globals.PARAM.angb1 * dRelSunhr)

    dTrans = dTotrad / dETRad
    dDfFrac = 1.0
    if dTrans >= 0.75:
        dDfFrac = 0.23
    elif dTrans < 0.75 and dTrans >= 0.35:
        dDfFrac = 1.33 - 1.46 * dTrans
    elif dTrans < 0.35 and dTrans >= 0.07:
        dA = dTrans - 0.07
        dDfFrac = 1.0 - 2.3 * dA * dA
   
    dfrad = dDfFrac * dTotrad
    drrad = dTotrad - dfrad
    return (dfrad,drrad)

# Hourly solar irradiance (W m-2). th (solar time),
#   dfrad (diffuse irradiance), drrad (direct irradiance)
#   where both dfrad and drrad will be set inside this function
def HourRad(th, dfrad, drrad):

    dfrad,drrad=DayRad(dfrad, drrad)
    dDayTotRad = dfrad + drrad
    dDecl = SolarDeclination()
    dA = math.sin(dDecl) * math.sin(globals.PARAM.lat)
    dB = math.cos(dDecl) * math.cos(globals.PARAM.lat)
    dAoB = dA / dB
    dPhi = (globals.Globals.gPI * dDayTotRad / 86400) / (dA * math.acos(-dAoB) + dB * math.sqrt(1 - dAoB * dAoB))
    dCoefA = -dB * dPhi
    dCoefB = dA * dPhi
    dTotrad = dCoefA * math.cos(globals.Globals.gPI * th / 12) + dCoefB
    dSinb = math.sin(SolarElevation(th))

    # determine fraction diffuse and direct:
    dDfFrac = 1.0 # fraction diffuse
    if dTotrad > 0.0 and dSinb >= 0.0:
        dR = 0.847 - 1.61 * dSinb + 1.04 * dSinb * dSinb
        dK = (1.47 - dR) / 1.66
        dHourET = HourETRad(th) # won't be zero
        dTrans = dTotrad / dHourET
        if dTrans > dK:
            dDfFrac = dR
        elif dTrans <= dK and dTrans > 0.35:
            dDfFrac = 1.47 - 1.66 * dTrans
        elif dTrans <= 0.35 and dTrans > 0.22:
            dA = dTrans - 0.22
            dDfFrac = 1.0 - 6.4 * dA * dA
    elif dTotrad < 0.0:
        dTotrad = 0.0

    dfrad = dDfFrac * dTotrad
    drrad = dTotrad - dfrad
    return (dfrad,drrad)

# Saturated vapor pressure (mbar). temp (temperature)
def Svp(temp):
    return (6.1078 * math.exp(17.269 * temp / (temp + 237.3)))

# Hourly vapor pressure (mbar). temp is air temperature (deg. C)
def HourVP(temp):
    # vapor pressure (convert RH to vapour pressure)
    return (wthrfile.DAILY_METEO.rh * Svp(temp) / 100.0) # assume constant RH

# Hourly wind speed (m s-1)
def HourWind():
    return wthrfile.DAILY_METEO.wind # assume constant wind speed

# Hourly temperature (deg. C). Simulation based on max and 
#   min daily air temperatures. th (time in hours)
def HourTemp(th):
    P = 1.5 # offset 1.5h after sunrise and noon
    dRise = Sunrise()
    dSet = Sunset()
    #print(dRise)
    #print(dSet)

    # how temperature is simulated depends on the current hour
    dTemp = 0.0
    if th >= (dRise+P) and th <= dSet:
        # diurnal
        dTau = globals.Globals.gPI * (th - dRise - P) / (dSet - dRise)
        dTemp = wthrfile.DAILY_METEO.tmin + (wthrfile.DAILY_METEO.tmax - wthrfile.DAILY_METEO.tmin) * math.sin(dTau)
    else:
        # nocturnal
        if th < (dRise+P):
            th = th + 24.0

        # temperature at time of sunset
        dTau = globals.Globals.gPI * (dSet - dRise - P) / (dSet - dRise)
        dTset = wthrfile.DAILY_METEO.tmin + (wthrfile.DAILY_METEO.tmax - wthrfile.DAILY_METEO.tmin)*math.sin(dTau)
        dSlope = (wthrfile.DAILY_METEO.tmin - dTset)/(dRise + P + 24 - dSet)
        dTemp = dTset + dSlope * (th - dSet)

    return dTemp

# Hourly net radiation. th is local solar time.
def HourNetRad(th):
    dDfrad = 0.0
    dDrrad = 0.0
    dDfrad,dDrrad=HourRad(th, dDfrad, dDrrad)
    dTotrad = 0.85 * (dDfrad + dDrrad) # net shortwave
    SB = 0.0000000567
    dTa = HourTemp(th) + 237.3
    dRLu = SB * dTa ** 4
    dRLd = 0.00000935 * SB * dTa ** 6
    dRLn = dRLd - dRLu # net longwave
    dRLn = dRLn * (0.2 + 0.8 * (wthrfile.DAILY_METEO.avsun / Daylength()))
    return (dTotrad + dRLn)


