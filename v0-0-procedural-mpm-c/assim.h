// header file: assim.h
#ifndef ASSIM_H
#define ASSIM_H

#include "et.h"

// Changes the given model parameter based on current temperature.
//   k25 (model parameter at 25 deg. C), q10 (Q10 value),
//   lftemp (leaf temperature, deg. C)
double Q10(double k25, double q10, double lftemp)
{
   return (k25 * pow(q10, (lftemp - 25.0) / 10.0));
}

// Rubisco capacity rate (umol m-2 s-1)
//   lftemp (leaf temperature, deg. C)
double Vcmax(double lftemp)
{
   double const VCMAX = 200.0;   // max rate
   double dExp = 1.0 + exp(0.128 * (lftemp - 40.0));
   return (Q10(VCMAX, 2.4, lftemp) / dExp);
}

// CO2 compensation point (umol mol-1)
//   lftemp (leaf temperature, deg. C)
double CO2CompensationPt(double lftemp)
{
   double const TAU = 2600.0;    // specificity factor
   double const OA = 210000.0;   // O2 concentration
   double dTau = Q10(TAU, 0.57, lftemp);
   return (OA * 0.5 / dTau);
}

// Gross photosynthesis limited by light (umol m-2 s-1)
//   par (PAR, umol m-2 s-1), lftemp (leaf temperature, deg. C),
//   ci (intercellular CO2 concentration, umol mol-1)
double LightLimited(double par, double lftemp, double ci)
{
   double const dQe = 0.06;   // quantum yield
   double const ALPHA = 0.8;  // absorption fraction
   double dT = CO2CompensationPt(lftemp);
   double dA = ALPHA * dQe * par * (ci - dT);
   double dB = ci + 2.0 * dT;
   return (dA / dB);
}

// Gross photosynthesis limited by Rubisco capacity (umol m-2 s-1)
//   lftemp (leaf temperature, deg. C),
//   ci (intercellular CO2 concentration, umol mol-1)
double RubiscoLimited(double lftemp, double ci)
{
   double const KCMAX = 300.0;         // O2
   double const KOMAX = 300000.0;      // CO2
   double const OA = 210000.0;         // O2 concentration
   double dKc = Q10(KCMAX, 2.1, lftemp);
   double dKo = Q10(KOMAX, 1.2, lftemp);
   double dKm = dKc * (1.0 + OA / dKo);
   double dA = Vcmax(lftemp) * (ci - CO2CompensationPt(lftemp));
   double dB = ci + dKm;
   return (dA / dB);
}

// Gross photosynthesis limited by sucrose sink (umol m-2 s-1)
//   lftemp (leaf temperature, deg. C)
double SinkLimited(double lftemp)
{
   return (0.5 * Vcmax(lftemp));
}

// Gross leaf photosynthesis (umol m-2 s-1)
//   par (PAR, umol m-2 s-1), lftemp (leaf temperature, deg. C),
double LeafAssim(double par, double lftemp)
{
   double const CI = 245.0; // internal CO2 concentration (umol mol-1)
   double dMin = RubiscoLimited(lftemp, CI);
   double dJc = LightLimited(par, lftemp, CI);
   double dJs = SinkLimited(lftemp);

   // find the most limiting factor to assimilation
   if (dMin > dJc)
      dMin = dJc;

   if (dMin > dJs)
      dMin = dJs;

   return dMin;
}

// Gross canopy photosynthesis (per unit ground area)
//   (umol CO2 m-2 s-1).
//   th (local solar time, hour), lai (leaf area index, m2 m-2),
//   lftemp (leaf temperature, deg. C)
double HourCanopyAssim(double th, double lftemp)
{
   // absorbed PAR
   double dQsl = 0.0, dQsh = 0.0;
   AbsorbedHourPAR(th, dQsl, dQsh);
   dQsl = 4.55 * dQsl;     // convert to umol m-2 s-1
   dQsh = 4.55 * dQsh;
   // leaf assimilation:
   double dAsl = LeafAssim(dQsl, lftemp);
   double dAsh = LeafAssim(dQsh, lftemp);
   // sunlit and shaded LAI
   double dLsl = 0.0, dLsh = 0.0;
   LAI(th, dLsl, dLsh);
   return (dAsl * dLsl + dAsh * dLsh); // in umol CO2 m-2 s-1
}

// Daily gross canopy photosynthesis (per unit ground area)
//   (umol CO2 m-2 day-1)
double DayCanopyAssim()
{
   // 5-point gauss integration over sunrise to sunset:
   double const ABS[5] = {0.0469101, 0.2307534, 0.5000000,
                          0.7692465, 0.9530899};
   double const WGT[5] = {0.1184635, 0.2393144, 0.2844444,
                          0.2393144, 0.1184635};

   double dTsr = Sunrise();
   double dIval = Sunset() - dTsr;
   double dAssim = 0.0;
   for (int i=0; i<5; ++i)
   {
      double dHour = dTsr + ABS[i] * dIval;   // current hour
      double dET = 0.0, dETs = 0.0, dETc = 0.0;
      double dH = 0.0, dHs = 0.0, dHc = 0.0;
      HourHeatFluxes(dHour, dET, dETs, dETc, dH, dHs, dHc);
      double dTf = CanopyTemp(dHour, dH, dHc);
      double dAn = HourCanopyAssim(dHour, dTf);
      // x 3600 to convert sec to hour:
      dAssim = dAssim + dAn * 3600.0 * WGT[i] * dIval;
   }

   return dAssim;
}

#endif
