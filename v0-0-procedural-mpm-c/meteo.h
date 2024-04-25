// header file: meteo.h
#ifndef METEO_H
#define METEO_H

#include <cmath>
#include "globals.h"

// Converts local time to local solar time. t (local time)
//double LocalSolarTime(double t)
//{
//   double dLonSM = (int)(12.0 * gParam.lon / gPI) * (gPI / 12);
//   double dLonCor = 12.0 * (dLonSM - gParam.lon) / gPI;
//   double dB = 2.0 * gPI * (gDoy - 81) / 364;
//   double dEoT = 9.87 * sin(2*dB) - 7.53 * cos(dB) - 1.5 * sin(dB);
//   double dTime = t + dLonCor + dEoT / 60;
//
//   // ensure time stays within a 24-hour period
//   if (dTime < 0.0)
//      dTime = 24.0 + dTime;
//   else if (dTime >= 24.0)
//      dTime = dTime - 24.0;
//
//   return dTime;
//}
//
//// Converts local solar time to local time. th (local solar time)
//double LocalTime(double th)
//{
//   double dLonSM = static_cast<int>((12.0 * gParam.lon / gPI))
//                   * (gPI / 12);
//   double dLonCor = 12.0 * (dLonSM - gParam.lon) / gPI;
//   double dB = 2.0 * gPI * (gDoy - 81) / 364;
//   double dEoT = 9.87 * sin(2 * dB) - 7.53 * cos(dB) - 1.5 * sin(dB);
//   double dTime = th - dLonCor - dEoT / 60;
//
//   // ensure time stays within a 24-hour period
//   if (dTime < 0.0)
//      dTime = 24.0 + dTime;
//   else if (dTime >= 24.0)
//      dTime = dTime - 24.0;
//
//   return dTime;
//}

// Solar declination (radians)
double SolarDeclination()
{
   return (-0.4093 * cos(2 * gPI * (gDoy + 10) / 365));
}

// Solar angle from horizontal (radians). Note: function will
//   return a -ve value if sun is below horizon. th (solar time)
double SolarElevation(double th)
{
   double dDecl = SolarDeclination();
   double dA = sin(dDecl) * sin(gParam.lat);
   double dB = cos(dDecl) * cos(gParam.lat);
   double dHa = gPI * (th - 12) / 12;
   return (asin(dA + dB * cos(dHa)));  // can be -ve
}

// Solar angle from north in an eastward direction (radians)
//   th (solar time)
double SolarAzimuth(double th)
{
   double dDecl = SolarDeclination();
   double dElev = SolarElevation(th);
   double dVal = (sin(gParam.lat) * sin(dElev) - sin(dDecl)) /
                 (cos(gParam.lat) * cos(dElev));

   if (dVal > 1.0)   // problem if dVal > 1 (see *below)
      dVal = 1.0;
  
   double dAzi = acos(dVal);  // * ensure dVal <= 1.0
   if (th < 12)
      dAzi = -dAzi;  // -ve for before solar noon, else +ve
   return (gPI + dAzi);
}

// Length of day (hours) - from sun rise to sun set
double Daylength()
{
   double dDecl = SolarDeclination();
   double dA = sin(dDecl) * sin(gParam.lat);
   double dB = cos(dDecl) * cos(gParam.lat);
   return (24 * acos(-dA / dB) / gPI);
}

// Solar time of sun set (hours)
double Sunset()
{
   double dDecl = SolarDeclination();
   double dA = sin(dDecl) * sin(gParam.lat);
   double dB = cos(dDecl) * cos(gParam.lat);
   return (12 * acos(-dA / dB) / gPI + 12);
}

// Solar time of sun rise (hours)
double Sunrise()
{
   return (24 - Sunset());
}

// Hourly extraterrestrial solar irradiance (W m-2). th (solar time)
double HourETRad(double th)
{
   double dE0 = 1 + 0.033 * cos(2 * gPI * (gDoy - 10) / 365);
   double dETRad = 1370 * dE0 * sin(SolarElevation(th));
   if (dETRad < 0.0)
      dETRad = 0.0;     // no -ve radiation values
   return dETRad;
}

// Daily extraterrestrial solar irradiance (J m-2 day-1)
double DayETRad()
{
   double dDecl = SolarDeclination();
   double dA = sin(dDecl) * sin(gParam.lat);
   double dB = cos(dDecl) * cos(gParam.lat);
   double dAoB = dA / dB;
   double dSinbe = (24.0 / gPI) * (dA * acos(-dAoB) +
                                   dB * sqrt(1 - dAoB * dAoB));
   double dE0 = 1.0 + 0.033 * cos(2 * gPI * (gDoy - 10) / 365);
   double dScc = 1370.0 * dE0;
   return (3600 * dScc * dSinbe);
}

// Daily solar irradiance (J m-2 day-1)
//   dfrad (diffuse irradiance), drrad (direct irradiance)
//   where both dfrad and drrad will be set inside this function
void DayRad(double &dfrad, double &drrad)
{
   double dRelSunhr = gMet.avsun / Daylength();
   double dETRad = DayETRad();
   double dTotrad = dETRad * (gParam.angb0 + gParam.angb1 * dRelSunhr);

   double dTrans = dTotrad / dETRad;
   double dDfFrac = 1.0;
   if (dTrans >= 0.75)
      dDfFrac = 0.23;
   else if (dTrans < 0.75 && dTrans >= 0.35)
      dDfFrac = 1.33 - 1.46 * dTrans;
   else if (dTrans < 0.35 && dTrans >= 0.07)
   {
      double dA = dTrans - 0.07;
      dDfFrac = 1.0 - 2.3 * dA * dA;
   }

   dfrad = dDfFrac * dTotrad;
   drrad = dTotrad - dfrad;
}

// Hourly solar irradiance (W m-2). th (solar time),
//   dfrad (diffuse irradiance), drrad (direct irradiance)
//   where both dfrad and drrad will be set inside this function
void HourRad(double th, double &dfrad, double &drrad)
{
   DayRad(dfrad, drrad);
   double dDayTotRad = dfrad + drrad;
   double dDecl = SolarDeclination();
   double dA = sin(dDecl) * sin(gParam.lat);
   double dB = cos(dDecl) * cos(gParam.lat);
   double dAoB = dA / dB;
   double dPhi = (gPI * dDayTotRad / 86400) /
                 (dA * acos(-dAoB) + dB * sqrt(1 - dAoB * dAoB));
   double dCoefA = -dB * dPhi;
   double dCoefB = dA * dPhi;
   double dTotrad = dCoefA * cos(gPI * th / 12) + dCoefB;
   double dSinb = sin(SolarElevation(th));

   // determine fraction diffuse and direct:
   double dDfFrac = 1.0;      // fraction diffuse
   if (dTotrad > 0.0 && dSinb >= 0.0)
   {
      double dR = 0.847 - 1.61 * dSinb + 1.04 * dSinb * dSinb;
      double dK = (1.47 - dR) / 1.66;
      double dHourET = HourETRad(th);     // won't be zero
      double dTrans = dTotrad / dHourET;
      if (dTrans > dK)
         dDfFrac = dR;
      else if (dTrans <= dK && dTrans > 0.35)
         dDfFrac = 1.47 - 1.66 * dTrans;
      else if (dTrans <= 0.35 && dTrans > 0.22)
      {
         double dA = dTrans - 0.22;
         dDfFrac = 1.0 - 6.4 * dA * dA;
      }
   }
   else if (dTotrad < 0.0)
      dTotrad = 0.0;

   dfrad = dDfFrac * dTotrad;
   drrad = dTotrad - dfrad;
}

// Saturated vapor pressure (mbar). temp (temperature)
double Svp(double temp)
{
   return (6.1078 * exp(17.269 * temp / (temp + 237.3)));
}

// Hourly vapor pressure (mbar). temp is air temperature (deg. C)
double HourVP(double temp)
{
   // vapor pressure (convert RH to vapour pressure)
   return (gMet.rh * Svp(temp) / 100.0); // assume constant RH
}

// Hourly wind speed (m s-1)
double HourWind()
{
   return gMet.wind;  // assume constant wind speed
}

// Hourly temperature (deg. C). Simulation based on max and 
//   min daily air temperatures. th (time in hours)
double HourTemp(double th)
{
   double const P = 1.5;   // offset 1.5h after sunrise and noon
   double dRise = Sunrise(), dSet = Sunset();

   // how temperature is simulated depends on the current hour
   double dTemp = 0.0;
   if (th >= (dRise+P) && th <= dSet)
   {
      // diurnal
      double dTau = gPI * (th - dRise - P) / (dSet - dRise);
      dTemp = gMet.tmin + (gMet.tmax - gMet.tmin) * sin(dTau);
   }
   else
   {
      // nocturnal
      if (th < (dRise+P))
         th = th + 24.0;

      // temperature at time of sunset
      double dTau = gPI * (dSet - dRise - P) / (dSet - dRise);
      double dTset = gMet.tmin + (gMet.tmax - gMet.tmin)*sin(dTau);
      double dSlope = (gMet.tmin - dTset)/(dRise + P + 24 - dSet);
      dTemp = dTset + dSlope * (th - dSet);
   }

   return dTemp;
}

// Hourly net radiation. th is local solar time.
double HourNetRad(double th)
{
   double dDfrad = 0.0, dDrrad = 0.0;
   HourRad(th, dDfrad, dDrrad);
   double dTotrad = 0.85 * (dDfrad + dDrrad);   // net shortwave
   double const SB = 0.0000000567;
   double dTa = HourTemp(th) + 237.3;
   double dRLu = SB * pow(dTa, 4);
   double dRLd = 0.00000935 * SB * pow(dTa, 6);
   double dRLn = dRLd - dRLu;    // net longwave
   dRLn = dRLn * (0.2 + 0.8 * (gMet.avsun / Daylength()));
   return (dTotrad + dRLn);
}

#endif
