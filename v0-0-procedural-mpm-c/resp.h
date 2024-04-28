// header file: resp.h
#ifndef RESP_H
#define RESP_H

#include "watbal.h"
#include "assim.h"

// Updates a given value (val) with its rate of change (rate)
void Intgrl(double &val, double rate)
{
   val = val + rate * gParam.delta;   // Euler's integration
}

// Increment to next development growth stage based on current
//   mean daily air temperature
void NextDVS()
{
   double dTemp = (gMet.tmin + gMet.tmax) / 2.0;   // mean
   double dDvr = 0.0;
   if (gParam.dvs <= 1.0)
      gParam.tabdvr01.LookUp(dTemp, dDvr); // pre-anthesis
   else
      gParam.tabdvr02.LookUp(dTemp, dDvr); // post-anthesis

   Intgrl(gParam.dvs, dDvr);
}

// Leaf area index (LAI) growth (m2 m-2).
void LAIGrowth()
{
   double dSla = 0.0;
   gParam.tabsla.LookUp(gParam.dvs, dSla);
   gParam.lai = gParam.dm.gnleaves * dSla;
}

// Increment rooting depth. limit is reduction to growth (unitless).
void RootDepthGrowth(double limit)
{
   // increment only if current depth has not reached maximum rooting
   //   depth, soil water content is less than wilting point, and
   //   only during pre-anthesis period
   if (gParam.len02 < gParam.maxlen &&
       gParam.vwc02 > gParam.vwcwp && gParam.dvs < 1.0)
      Intgrl(gParam.len02, limit * gParam.rootgrate);
}

// Increment plant height growth.
//   limit is reduction to growth (unitless).
void PlantHeightGrowth(double limit)
{
   // use different base temperature if for post-anthesis
   double dBaseTemp = gParam.basetemp01;
   if (gParam.dvs > 1.0)
      dBaseTemp = gParam.basetemp02;

   double dTemp = (gMet.tmin + gMet.tmax) / 2.0;   // mean
   double dTs = dTemp - dBaseTemp;   // temperature sum for the day
   if (dTs < 0.0)
      dTs = 0.0;        // no temp. sum if less than base temp.

   Intgrl(gTs, dTs);    // increment total temperature sum
   double dA = dTs * gParam.hgtb1 * gParam.hgtb0 *
               gParam.hgtmax * exp(-gParam.hgtb1 * gTs);
   double dB = 1.0 + gParam.hgtb0 * exp(-gParam.hgtb1 * gTs);
   double dRate = limit * dA / (dB * dB);    // daily height growth
   Intgrl(gParam.hgt, dRate);
}

// Maintenance respiration: determines the assimilates used for
//   maintenance (g CH20 m-2 day-1).
void RespMaint()
{
   // respiration coefficients (g CH2O g-1 DM day-1):
   double const MAINTST = 0.015, MAINTLV = 0.030,
                MAINTRT = 0.015, MAINTSO = 0.010;
   double dMaint = MAINTST * gParam.dm.stem +
                   MAINTLV * gParam.dm.gnleaves +
                   MAINTRT * gParam.dm.roots +
                   MAINTSO * gParam.dm.storage;
   double dTemp = (gMet.tmax + gMet.tmin) / 2.0;   // mean
   dMaint = Q10(dMaint, 2.0, dTemp);   // temperature correction
   // need to correct for decrease in metabolic activity with age:
   double dTotal = gParam.dm.gnleaves + gParam.dm.ddleaves;
   dMaint = dMaint * gParam.dm.gnleaves / dTotal;
   if (gAssim.gphot < dMaint)
      dMaint = gAssim.gphot; // not enough assimilates for maintenance

   gAssim.maint = dMaint;    // store it
}

// Determines daily death rate for leaves (day-1)
double LeafDeathRate()
{
   double dDDAge = 0.0;       // leaf death due to age
   if (gParam.dvs > 1.0)      // death only after anthesis
   {  
      double dN = 2.0 - gParam.dvs;
      if (dN < 0.1)
         dN = 0.1;

      double dDvr = 0.0;
      double dTemp = (gMet.tmin + gMet.tmax) / 2.0;
      gParam.tabdvr02.LookUp(dTemp, dDvr);
      dDDAge = dDvr / dN;
   }

   double dDDShade = 0.0;     // leaf death due to self-shading
   double const LAICR = 4.0;  // critical LAI for self-shading death
   if (gParam.lai > LAICR)
   {
      double const DDLEAF = 0.03;   // max. dead leaf coefficient
      dDDShade = DDLEAF * (gParam.lai - LAICR) / LAICR;
      if (dDDShade > DDLEAF)
         dDDShade = DDLEAF;
   }

   double dRdr = dDDAge;
   if (dDDAge < dDDShade)  // select the maximum death rate
      dRdr = dDDShade;
   
   return dRdr;
}

// Growth respiration: determines the assimilates used for
//   growth (g CH20 m-2 day-1).
void RespGrowth()
{
   gAssim.growth = gAssim.gphot - gAssim.maint; // leftover for growth

   // read in tabulated data for DM partitioning:
   double dPartst = 0.0, dPartlv = 0.0, dPartrt = 0.0, dPartso = 0.0;
   gParam.tabst.LookUp(gParam.dvs, dPartst);
   gParam.tablv.LookUp(gParam.dvs, dPartlv);
   gParam.tabrt.LookUp(gParam.dvs, dPartrt);
   gParam.tabso.LookUp(gParam.dvs, dPartso);

   // growth respiration coefficients (g CH20 g-1 DM):
   double const GROST = 1.513, GROLV = 1.463,
                GRORT = 1.444, GROSO = 1.415;
   double dGroTot = GROST * dPartst + GROLV * dPartlv +
                    GRORT * dPartrt + GROSO * dPartso;
   if (dGroTot > 0.0)
      dGroTot = gAssim.growth / dGroTot;  // in g DM m-2 day-1

   double dGst = dPartst * dGroTot;       // stem
   double dGlv = dPartlv * dGroTot;       // leaves
   double dGrt = dPartrt * dGroTot;       // roots
   double dGso = dPartso * dGroTot;       // storage organs
   double dGdlv = gParam.dm.gnleaves * LeafDeathRate();  // leaf death

   // update all DM weights:
   Intgrl(gParam.dm.gnleaves, dGlv - dGdlv);
   Intgrl(gParam.dm.ddleaves, dGdlv);
   Intgrl(gParam.dm.stem, dGst);
   Intgrl(gParam.dm.roots, dGrt);
   Intgrl(gParam.dm.storage, dGso);
}

// Determines the assimilates produced and used for current
//   development stage. Returns the previous parameter (PARAM) values.
//   potT  - potential transpiration (mm day-1).
void Grow(double potT)
{
   double dLimit = GrowthReduction(potT); // water stress

   // assimilates produced:
   gAssim.gphot = 0.0;
   double dGphot = DayCanopyAssim();
   dGphot = dGphot * dLimit;  // reduced growth due to water stress
   if (dGphot > 0.0)
      gAssim.gphot = dGphot * 30.0 / 1000000;  // in g CH2O m-2 day-1

   // assimilates used:
   RespMaint();
   RespGrowth();

   LAIGrowth();                // update LAI
   PlantHeightGrowth(dLimit);  // update plant height
   RootDepthGrowth(dLimit);    // update rooting depth
}

#endif
