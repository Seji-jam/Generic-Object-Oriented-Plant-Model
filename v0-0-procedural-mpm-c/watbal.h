// header file: watbal.h
#ifndef WATBAL_H
#define WATBAL_H

#include <cmath>
#include "globals.h"

// Water content before redistribution (mm)
//   length = soil layer thickness (m)
//   wc = current water content (m3 m-3)
//   wcperc = water percolating from the above soil layer (mm day-1)
double WcBeforeRedist(double length, double wc, double wcperc)
{
   double dC = 1000.0 * length;  // for conversion
   double dWcsat = gParam.vwcsat * dC;
   double dPercEx = 0.0;         // percolation to below
   double dStore = wc + wcperc;  // total water amount
   if (dStore > dWcsat)
      dPercEx = dStore - dWcsat; // total exceeds saturation
   return(wc + wcperc - dPercEx);
}

// Water content after redistribution (mm)
//   length = soil layer thickness (m)
//   wc = current water content (m3 m-3)
double WcAfterRedist(double length, double wc)
{
   double dC = 1000.0 * length;  // for conversion
   double dVwc = wc / dC;        // want unit in m3 m-3
   double dExp = exp((gParam.alpha / gParam.vwcsat) *
                     (gParam.vwcsat - dVwc));
   double dVal = (gParam.alpha * gParam.ksat * gParam.delta) /
                 (length * gParam.vwcsat);
   double dLog = log(dVal + dExp);
   double dWc = gParam.vwcsat - (gParam.vwcsat/gParam.alpha) * dLog;
   if (dWc < 0)
      dWc = 0.0;        // no negative water content
   return (dWc * dC);   // convert back to mm
}

// Actual soil evaporation (mm day-1)
//   potE = potential soil evaporation (mm day-1)
//   e01 and e02 = actual evaporation from layer 1 and 2 (mm day-1)
//                 (set within this function)
void ActualE(double potE, double &e01, double &e02)
{
   double dPow = pow((3.6073 * gParam.vwc01 / gParam.vwcsat),
                     -9.3172);
   double dReduction = 1.0 / (1.0 + dPow);
   e01 = potE * dReduction * 0.26;  // 26% from top layer
   e02 = potE * dReduction * 0.74;  // rest from root zone
}

// Actual plant transpiration (mm day-1)
//   p = 0.5 or 0.3 for C3 and C4 plants, respectively
//   potT = potential soil transpiration (mm day-1)
//   t01 and t02 = actual transpiration from layer 1 and 2 (mm day-1)
//                 (set within this function)
void ActualT(double p, double potT, double &t01, double &t02)
{
   // critical water content
   double dVwccr = gParam.vwcwp + 
                   p * (gParam.vwcsat - gParam.vwcwp);
   double dReduction = (gParam.vwc02 - gParam.vwcwp) /
                       (dVwccr - gParam.vwcwp);
   if (dReduction > 1.0)
      dReduction = 1.0;

   t01 = 0.0;                 // no active roots in top layer
   t02 = potT * dReduction;   // all roots in root zone
}

// Final volumetric water content (m3 m-3)
//   length = soil layer thickness (m)
//   vwc = current volumetric water content (m3 m-3)
//   ea and ta = actual evaporation and transpiration (mm day-1)
//   wcperc0 and wcperc1 = percolation from above and to below
//                         (mm day-1)
//   Note: wcperc1 is set within this function.
double Vwc(double length, double vwc, double ea, double ta,
           double wcperc0, double &wcperc1)
{
   double dC = 1000.0 * length;     // for conversions
   double dWc = vwc * dC;           // want unit in mm
   dWc = WcBeforeRedist(length, dWc, wcperc0);
   double dBal = WcAfterRedist(length, dWc);
   dBal = dBal - ea - ta;           // final water amount
   if (dBal < 0.0)
      dBal = 0.0;

   wcperc1 = dWc - ea - ta - dBal;  // percolation to below
   if (wcperc1 < 0.0)
      wcperc1 = 0.0;

   return (dBal / dC);              // unit m3 m-3
}

// Final volumetric water content (m3 m-3) for the two soil
//   layers. potE and potT are the potential evaporation and
//   transpiration (mm day-1).
void DailyWaterContent(double potE, double potT)
{
   double ea01 = 0.0, ea02 = 0.0;
   ActualE(potE, ea01, ea02);
   double ta01 = 0.0, ta02 = 0.0;
   ActualT(0.5, potT, ta01, ta02);  // assume a C3 plant
   double dPerc01 = 0.0;
   gParam.vwc01 = Vwc(gParam.len01, gParam.vwc01, ea01, ta01,
                      gMet.rain, dPerc01);
   gParam.vwc02 = Vwc(gParam.len02, gParam.vwc02, ea02, ta02,
                      dPerc01, dPerc01);
}

// Fraction of growth reduced due to limited water. potT is potential
//   transpiration (mm day-1).
double GrowthReduction(double potT)
{
   double ta01 = 0.0, ta02 = 0.0;
   ActualT(0.5, potT, ta01, ta02);  // C3 plant
   return (ta02 / potT);            // fraction (0 to 1)
}

#endif
