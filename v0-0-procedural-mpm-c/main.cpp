// file: main.cpp
#include <iostream>
#include "resp.h"

int main()
{
   try
   {
      ifstream in("input.txt");   // input/data file (text)
      ofstream out("output.txt"); // output/results file

      in >> gParam;  // read model parameters from data file

      // set simulation clock
      clck now(gParam.year, gParam.month, gParam.day);

      // Simulation loop. Ends when either: 1) max. development stage
      //   (gParam.dvsstop) is exceeded, or 2) no. of simulation
      //   steps (gParam.steps) is reached -- whichever comes first.
      bool bLastRun = false;
      for (int i = 0; i <= gParam.steps && !bLastRun;
           ++i, now.Advance(gParam.delta), NextDVS())
      {
         // if max. development stage is exceeded, do the simulation
         //   for one last time to obtain and print last set of
         //   values, then end simulation (bLastRun set to true).
         bLastRun = (gParam.dvs >= gParam.dvsstop);

         wthrfile wf(gParam.fname, now);  // weather file
         ReadWeather(wf);                 // read weather values

         // output current parameter values:
         out << setw(10) << (i*gParam.delta) << "\t"   // elapsed days
             << setw(10) << gParam.dvs << "\t"
             << setw(10) << gParam.lai << "\t"
             << setw(10) << gParam.dm.gnleaves << "\t"
             << setw(10) << gParam.dm.ddleaves << "\t"
             << setw(10) << gParam.dm.stem << "\t"
             << setw(10) << gParam.dm.roots << "\t"
             << setw(10) << gParam.dm.storage << "\t"
             << setw(10) << gParam.hgt << "\t"
             << setw(10) << gParam.vwc01 << "\t"
             << setw(10) << gParam.vwc02 << "\t"
             << setw(10) << gParam.len02 << "\t";

         // daily heat fluxes:
         double dET = 0.0, dETs = 0.0, dETc = 0.0;
         double dH = 0.0, dHs = 0.0, dHc = 0.0;
         DayHeatFluxes(dET, dETs, dETc, dH, dHs, dHc);

         // daily soil water content (in mm day-1):
         double dPotE = dETs * 1000.0 / (2454000.0 * 998.0);
         double dPotT = dETc * 1000.0 / (2454000.0 * 998.0);
         DailyWaterContent(dPotE, dPotT);

         // obtain actual evapotranspiration:
         double dEa01 = 0.0, dEa02 = 0.0;
         ActualE(dPotE, dEa01, dEa02);
         double dTa01 = 0.0, dTa02 = 0.0;
         ActualT(0.5, dPotT, dTa01, dTa02);

         Grow(dPotT);         // plant growth (next stage)

         // output simulation results:
         out << setw(10) << gAssim.gphot << "\t"
             << setw(10) << gAssim.maint << "\t"
             << setw(10) << gAssim.growth << "\t"
             << setw(10) << (dETs/1000000) << "\t"  // in MJ m-2 day-1
             << setw(10) << (dETc/1000000) << "\t"
             << setw(10) << ((dETs + dETc)/1000000) << "\t"
             << setw(10) << (dHs/1000000) << "\t"
             << setw(10) << (dHc/1000000) << "\t"
             << setw(10) << ((dHs + dHc)/1000000) << "\t"
             << setw(10) << dPotE << "\t"
             << setw(10) << dPotT << "\t"
             << setw(10) << (dPotE + dPotT) << "\t"
             << setw(10) << (dEa01 + dEa02) << "\t"
             << setw(10) << (dTa01 + dTa02) << "\t"
             << setw(10) << (dEa01 + dEa02 + dTa01 + dTa02) << "\n";
      }
   }
   catch (exception &e)
   {
      cout << "Error occurred. " << e.what() << endl;
   }

   cout << "Simulation ended." << endl;
   return 0;
}
