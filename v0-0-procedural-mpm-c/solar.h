// header file: solar.h
#ifndef SOLAR_H
#define SOLAR_H

#include "meteo.h"

// Extinction coefficient for direct fluxes. th is solar time
double Kdr(double th)
{
   double dElev = SolarElevation(th);
   if (dElev < 0.00000001)
      return 0.0;    // sun is below horizon
   return (0.5 / sin(dElev));
}

// Extinction coefficient for diffuse fluxes
double Kdf()
{
   double dL = sqrt(gParam.lai);
   return ((1.0 + 0.1174 * dL) / (1.0 + 0.3732 * dL));
}

// Hourly intercepted of total solar radiation (W m-2)
//   th (solar time),
//   dfrad (diffuse irradiance), drrad (direct irradiance)
//   where both dfrad and drrad will be set inside this function
void InterceptHourRad(double th, double &dfrad, double &drrad)
{
   double const P = 0.11;        // reflection coefficient
   double const S = sqrt(0.50);  // scatter correction
   double dDfRad = 0.0, dDrRad = 0.0;
   HourRad(th, dDfRad, dDrRad);
   drrad = (1.0 - P) * dDrRad * (1.0 - exp(-Kdr(th) * S * gParam.lai));
   dfrad = (1.0 - P) * dDfRad * (1.0 - exp(-Kdf() * S * gParam.lai));
}

// Daily intercepted of total solar radiation (J m-2 day-1)
//   dfrad (diffuse irradiance), drrad (direct irradiance)
//   where both dfrad and drrad will be set inside this function
void InterceptDayRad(double &dfrad, double &drrad)
{
   // numerical integration using 5-point Gaussian method
   double const ABS[5] = {0.0469101, 0.2307534, 0.5000000,
                          0.7692465, 0.9530899};
   double const WGT[5] = {0.1184635, 0.2393144, 0.2844444,
                          0.2393144, 0.1184635};
   double dTsr = Sunrise();
   double dIval = Sunset() - dTsr;
   dfrad = 0.0; drrad = 0.0;
   for (int i=0; i<5; ++i)
   {
      double dHr = dTsr + ABS[i] * dIval;
      double dDfrad = 0.0, dDrrad = 0.0;
      InterceptHourRad(dHr, dDfrad, dDrrad);
      // x 3600 to convert sec to hour:
      dfrad = dfrad + dDfrad * 3600 * WGT[i] * dIval;
      drrad = drrad + dDrrad * 3600 * WGT[i] * dIval;
   }
}

// Absorbed PAR (W m-2) by sunlit and shaded leaves.
//   th (solar time), 
//   both sunlit and shaded will be set inside this function
void AbsorbedHourPAR(double th, double &sunlit, double &shaded)
{
   double const P = 0.04;        // reflection coefficient
   double const A = 0.80;        // scatter coefficient
   double const S = sqrt(A);     // scatter correction
   double dDfPar = 0.0, dDrPar = 0.0;
   HourRad(th, dDfPar, dDrPar);
   
   dDfPar = 0.5 * dDfPar;        // 50% total radiation is PAR
   dDrPar = 0.5 * dDrPar;

   double dKdr = Kdr(th);
   double dI = (1.0 - P) * dDrPar;
   double dIpdr = dI * exp(-dKdr * S * gParam.lai);// total direct
   double dIpdrdr = dI * exp(-dKdr * gParam.lai);  // direct of direct
   double dIpdra = (dIpdr - dIpdrdr) / 2;          // scatter beams
   double dN = Kdf() * S * gParam.lai;
   dI = (1.0 - P) * dDfPar;
   double dIpdf = dI * (1.0 - exp(-dN)) / dN;   // diffuse

   sunlit = A * (dKdr * dDrPar + dIpdf + dIpdra);
   shaded = A * (dIpdf + dIpdra);
}

// Leaf area index for sunlit and shaded leaves.
//   th (solar time).
//   Both sunlit and shaded will be set inside this function
void LAI(double th, double &sunlit, double &shaded)
{
   sunlit = 0.0;
   double dKdr = Kdr(th);
   if (dKdr > 0.0)
      sunlit = (1.0 - exp(-dKdr * gParam.lai)) / dKdr;
   shaded = gParam.lai - sunlit;
}

#endif
