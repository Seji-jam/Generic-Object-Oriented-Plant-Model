// header file: et.h
#ifndef ET_H
#define ET_H

#include "solar.h"

// Slope of saturated vapor pressure with temperature (mbar/K).
//   temp (temperature)
double SlopeSvp(double temp)
{
   double dA = 25029.4 * exp(17.269 * temp / (temp + 237.3));
   double dB = pow(temp + 237.3, 2);
   return (dA / dB);
}

// Vapor pressure deficit at reference height (mbar).
//   temp (temperature)
double VpdRef(double temp)
{
   return (Svp(temp) - HourVP(temp));
}

// Vapor pressure deficit at mean canopy flow (mbar).
//   a (total energy available), et (total ET),
//   raa (aerodynamic resistance between mean canopy flow
//   and reference height), temp (temperature)
double VpdMcf(double a, double et, double raa, double temp)
{
   double dDelta = SlopeSvp(temp);
   double dA = (raa / 1221.09) *
               (dDelta * a - (dDelta + 0.658) * et);
   return (VpdRef(temp) + dA);
}

// Sensible heat for soil (W m-2).
//   rsa (aerodynamic resistance between mean canopy flow
//   and soil), rss (soil resistance), vpdmcf (vapor pressure
//   deficit at mean canopy flow), as (energy available to soil),
//   slopesvp (slope of saturated vapor pressure
//   with temperature)
double Hs(double rsa, double rss, double vpdmcf, double as,
          double slopesvp)
{
   double dA = 0.658 * as * (rss + rsa) - 1221.09 * vpdmcf;
   double dB = slopesvp * rsa + 0.658 * (rss + rsa);
   return (dA / dB);
}

// Sensible heat for crop (W m-2).
//   rca (aerodynamic resistance between mean canopy flow
//   and canopy), rcs (canopy resistance), vpdmcf (vapor pressure
//   deficit at mean canopy flow), ac (energy available to
//   canopy), slopesvp (slope of saturated vapor pressure
//   with temperature)
double Hc(double rca, double rcs, double vpdmcf, double ac,
          double slopesvp)
{
   double dA = 0.658 * ac * (rcs + rca) - 1221.09 * vpdmcf;
   double dB = slopesvp * rca + 0.658 * (rcs + rca);
   return (dA / dB);
}

// Latent heat for soil (W m-2).
//   rsa (aerodynamic resistance between mean canopy flow
//   and soil), rss (soil resistance), vpdmcf (vapor pressure
//   deficit at mean canopy flow), as (energy available to soil),
//   slopesvp (slope of saturated vapor pressure
//   with temperature)
double ETs(double rsa, double rss, double vpdmcf, double as,
           double slopesvp)
{
   double dA = slopesvp * as + (1221.09 * vpdmcf / rsa);
   double dB = slopesvp + 0.658 * (rss + rsa) / rsa;
   return (dA / dB);
}

// Latent heat for crop (W m-2).
//   rca (aerodynamic resistance between mean canopy flow
//   and canopy), rcs (canopy resistance), vpdmcf (vapor pressure
//   deficit at mean canopy flow), ac (energy available to
//   canopy), slopesvp (slope of saturated vapor pressure
//   with temperature)
double ETc(double rca, double rcs, double vpdmcf, double ac,
           double slopesvp)
{
   double dA = slopesvp * ac + (1221.09 * vpdmcf / rca);
   double dB = slopesvp + 0.658 * (rcs + rca) / rca;
   return (dA / dB);
}

// Soil heat flux (W m-2).
//   th (local solar time), rsn (solar radiation reaching the soil)
double G(double th, double rsn)
{
   double dSolarInc = (gPI/2.0) - SolarElevation(th);
   return (0.35 * cos(dSolarInc) * rsn);
}

// Zeroplane displacement (m)
double Zeroplane()
{
   return (0.64 * gParam.hgt);
}

// Surface roughness length (m)
double RoughLen()
{
   return (0.13 * gParam.hgt);
}

// Friction velocity (m2/s2). hgt (crop height),
//   windspd (wind speed at canopy top)
double FrictionVelocity(double windspd)
{
   double dD = Zeroplane();
   double dZ = RoughLen();
   return (0.4 * windspd / log((gParam.hgt - dD) / dZ));
}

// Calculate the two aerodynamic resistances (s/m):
//   rsa (aerodynamic resistance between mean canopy flow
//   and soil), and raa (aerodynamic resistance between mean
//   canopy flow and reference height), 
//   windspd (wind speed at canopy top)
void AeroResistances(double windspd, double &rsa, double &raa)
{
   double dD = Zeroplane();
   double dZ = RoughLen();
   double dFricVel = FrictionVelocity(windspd);
   double dKh = 0.4 * dFricVel * gParam.hgt;
   double dP = exp(-gParam.keddy * 0.004 / gParam.hgt) - 
               exp(-gParam.keddy * (dZ + dD) / gParam.hgt);
   rsa = (gParam.hgt * exp(gParam.keddy) / (gParam.keddy * dKh)) * dP;
   dP = exp(gParam.keddy * (1.0 - (dZ + dD) / gParam.hgt)) - 1.0;
   dP = (gParam.hgt / (gParam.keddy * dKh)) * dP;
   raa = log((gParam.refhgt - dD) / (gParam.hgt - dD)) /
         (0.4 * dFricVel) + dP;
}

// Boundary layer resistance (s/m), rca (resistance between
//   mean canopy flow and canopy).
//   windspd (wind speed at canopy top)
double BoundLayerResistance(double windspd)
{
   double dD = Zeroplane();
   double dZ = RoughLen();
   double dFricVel = FrictionVelocity(windspd);
   double dUh = (dFricVel / 0.4) * log ((gParam.hgt - dD) / dZ);
   double dRca = 0.012 * gParam.lai * (1.0 - exp(-gParam.kwind/2.0));
   dRca = gParam.kwind / (dRca * sqrt(dUh / gParam.leafwidth));
   return dRca;
}

// Canopy resistance (s/m), rcs.
//   totrad (total solar irradiance)
double CanopyResistance(double totrad)
{
   if (totrad <= 0.0)
      totrad = 0.01;   // avoid dividing by zero
  
   double dPar = 0.5 * totrad;   // PAR is 50% of solar rad.
   double dRcs = (gParam.rstA1 + dPar) / (gParam.rstA2 * dPar);

   double const LAICR = 0.5 * 4.0;  // critical LAI
   double dLe = gParam.lai;
   if (gParam.lai > LAICR)
      dLe = LAICR;

   return (dRcs / dLe);
}

// Soil resistance (s/m), rss
double SoilResistance()
{
   double dRssdry = (2.0 * 0.02) / (gParam.porosity * 0.0000247);
   double dExp = exp(-(1.0 / gParam.poredist) *
                      (gParam.vwc01 / gParam.vwcsat));
   return (dRssdry * dExp);
}

// Net energy available to soil and canopy (W m-2).
//   th (local solar time), lai (leaf area index),
//   as (energy available to soil), ac (energy available
//   to canopy). Both as and ac will be set inside this function.
void EnergySupply(double th, double &as, double &ac)
{
   ac = as = 0.0;
   double dRn = HourNetRad(th);
   if (dRn > 0.0)
   {
      double const S = sqrt(0.50);  // scatter correction
      double dRsn = dRn * exp(-Kdr(th) * S * gParam.lai);  // at soil
      ac = dRn - dRsn;
      as = dRn - G(th, dRsn) - ac;
   }
}

// Canopy temperature (deg. C).
//   th (local solar time),
//   h and hc - sensible heat for total and canopy
double CanopyTemp(double th, double h, double hc)
{
   double dWind = HourWind();
   double dRsa = 0.0, dRaa = 0.0;
   AeroResistances(dWind, dRsa, dRaa);
   double dT0 = ((h / 1221.09) * dRaa) + HourTemp(th);
   double dRca = BoundLayerResistance(dWind);
   double dTf = ((hc / 1221.09) * dRca) + dT0;
   return dTf;
}

// Latent and sensible heat fluxes (W m-2) at the given hour.
//   th (local solar time).
//   Variables set inside this function:
//   et, ets and etc - latent heat for total, soil and canopy
//   h, hs and hc - sensible heat for total, soil and canopy
void HourHeatFluxes(double th,double &et, double &ets, double &etc,
                    double &h, double &hs, double &hc)
{
   // hourly meteorological properties: /////////////////////////
   double dWind = HourWind();          // wind speed
   double dTemp = HourTemp(th);        // air temperature
   double dDfrad = 0.0, dDrrad = 0.0;
   HourRad(th, dDfrad, dDrrad);
   double dTotrad = dDfrad + dDrrad;   // solar irradiance

   // resistances: //////////////////////////////////////////////
   double dRsa = 0.0, dRaa = 0.0;
   AeroResistances(dWind, dRsa, dRaa);
   double dRca = BoundLayerResistance(dWind);
   double dRcs = CanopyResistance(dTotrad);
   double dRss = SoilResistance();

   // energy available to soil and plant: ///////////////////////
   double dAs = 0.0, dAc = 0.0;
   EnergySupply(th, dAs, dAc);

   // calculate fluxes now: /////////////////////////////////////
   double dDelta = SlopeSvp(dTemp);
   double dRa = (dDelta + 0.658) * dRaa;
   double dRc = (dDelta + 0.658) * dRca + 0.658 * dRcs;
   double dRs = (dDelta + 0.658) * dRsa + 0.658 * dRss;
   double dCc = 1.0 + (dRc * dRa / (dRs * (dRc + dRa)));
   dCc = 1.0 / dCc;
   double dCs = 1.0 + (dRs * dRa / (dRc * (dRs + dRa)));
   dCs = 1.0 / dCs;
   double dA = dAs + dAc;
   double dD = VpdRef(dTemp);

   double dP = dDelta * dA + 
               ((1221.09*dD - dDelta*dRca*dAs) / (dRaa + dRca));
   double dQ = dDelta + 0.658 * (1.0 + dRcs / (dRaa + dRca));
   double dPMc = dP / dQ;

   dP = dDelta * dA + 
        ((1221.09*dD - dDelta*dRsa*dAc) / (dRaa + dRsa));
   dQ = dDelta + 0.658 * (1.0 + dRss / (dRaa + dRsa));
   double dPMs = dP / dQ;

   et = dCc * dPMc + dCs * dPMs;
   double dD0 = VpdMcf(dA, et, dRaa, dTemp);
   ets = ETs(dRsa, dRss, dD0, dAs, dDelta);
   etc = et - ets;

   hs = Hs(dRsa, dRss, dD0, dAs, dDelta);
   hc = Hc(dRca, dRcs, dD0, dAc, dDelta);
   h = hs + hc;
}

// Latent and sensible heat fluxes (J m-2 day-1) for the day.
//   Variables set inside this function:
//   et, ets and etc - latent heat for total, soil and canopy
//   h, hs and hc - sensible heat for total, soil and canopy
void DayHeatFluxes(double &et, double &ets, double &etc,
                   double &h, double &hs, double &hc)
{
   // 5-point gauss integration:
   double const ABS[5] = {0.0469101, 0.2307534, 0.5000000,
                          0.7692465, 0.9530899};
   double const WGT[5] = {0.1184635, 0.2393144, 0.2844444,
                          0.2393144, 0.1184635};

   double dIval = 24.0;    // integration over 24 hours
   et = ets = etc = h = hs = hc = 0.0;       // reset to zeros
   for (int i=0; i<5; ++i)
   {
      double dHour = 0.0 + ABS[i] * dIval;   // current hour
      double det = 0.0, dets = 0.0, detc = 0.0;
      double dh = 0.0, dhs = 0.0, dhc = 0.0;
      HourHeatFluxes(dHour, det, dets, detc, dh, dhs, dhc);

      // x 3600 to convert sec to hour:
      et = et + det * 3600.0 * WGT[i] * dIval;
      ets = ets + dets * 3600.0 * WGT[i] * dIval;
      etc = etc + detc * 3600.0 * WGT[i] * dIval;
      h = h + dh * 3600.0 * WGT[i] * dIval;
      hs = hs + dhs * 3600.0 * WGT[i] * dIval;
      hc = hc + dhc * 3600.0 * WGT[i] * dIval;
   }
}

#endif
