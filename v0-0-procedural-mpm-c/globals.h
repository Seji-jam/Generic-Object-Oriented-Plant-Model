// header file: globals.h
#ifndef GLOBALS_H
#define GLOBALS_H

#include "wthrfile.h"   // weather file
#include "ini.h"        // initialization file

using namespace std;

// dry matter weight (all in g m-2)
struct WGT
{
   double stem, gnleaves, ddleaves, roots, storage;
};

// assimilates (all in g CH2O m-2 day-1)
struct ASSIM
{
   double gphot;     // gross photosynthesis
   double maint;     // maintenance respiration
   double growth;    // growth respiration
};

// the model parameters
struct PARAM
{
   // simulation and site properties:
   int delta;             // days per time step (days)
   int steps;             // no. of time steps to stop
   double dvsstop;        // growth development stage to stop
   int year, month, day;  // simulation date
   string fname;          // template for weather file name
   double lat, lon;       // site lat. and longitude (radians)
   double angb0, angb1;   // Angstrom coefficients

   // flux properties:
   double refhgt;         // reference height (m)
   double kwind;          // attenuation coefficient wind speed
   double keddy;          // attenuation coefficient eddy diffusivity
   double rstA1, rstA2;   // stomatal resist. coefficients
   double leafwidth;      // average leaf width (m)

   // soil properties:
   double poredist;       // pore size distribution index
   double porosity;       // porosity (fraction)
   double vwcwp;          // wilting point (m3 m-3)
   double vwcsat;         // saturation (m3 m-3)
   double vwc01, vwc02;   // water content (m3 m-3)
   double len01, len02;   // layer thickness (m)
   double maxlen;         // max. rooting depth (m)
   double alpha;          // slope of hydraulic conductivity
   double ksat;           // saturated hydraulic conductivity (m s-1)

   // plant properties:
   double dvs;            // growth development stage
   double basetemp01;     // pre-anthesis base growth temp. (deg. C)
   double basetemp02;     // post-anthesis base growth temp. (deg. C)
   double lai;            // leaf area index LAI (m2 m-2)
   double hgt;            // plant height (m)
   double hgtmax;         // max. plant height (m)
   double hgtb0;          // growth of plant height intercept
   double hgtb1;          // growth of plant height slope
   double rootgrate;      // roots growing rate (mm day-1)
   WGT dm;                // dry matter weights (g m-2)

   // tabulated data:
   afgen tabdvr01;   // temperature-development stage (pre-anthesis)
   afgen tabdvr02;   // temperature-development stage (post-anthesis)
   afgen tabst;      // DM partitioning to stem
   afgen tablv;      // DM partitioning to leaves
   afgen tabrt;      // DM partitioning to roots
   afgen tabso;      // DM partitioning to storage organs
   afgen tabsla;     // specific leaf area
};

// global variables:
double const gPI = 3.1415926536;
PARAM gParam;            // model parameters
DAILY_METEO gMet;        // daily weather properties
ASSIM gAssim;            // daily assimilates produced and used
int gDoy = 0;            // current day of year (e.g., Jan 1 = 1)
double gTs = 0.0;        // temperature (heat) sum

// Reads weather file and sets the meteorological parameters.
//    Also sets the day of year (gDoy) and date variables
//    to current date.
void ReadWeather(wthrfile &wf)
{
   if (wf.SetDailyMeteo())
   {
      // weather file read successfully: set the global variables
      gMet = wf.DailyMeteo();
      gDoy = wf.Now().DayOfYear();
      gParam.year = wf.Now().Date().year;
      gParam.month = wf.Now().Date().month;
      gParam.day = wf.Now().Date().day;
   }
}

// Sets the model parameters from data file
ifstream &operator>>(ifstream &in, PARAM &prm)
{
   ini init(in);

   // simulation and site properties:
   init.Set("Interval", prm.delta);
   init.Set("TimeStepsStop", prm.steps);
   init.Set("DevStageStop", prm.dvsstop);
   init.Set("Year", prm.year);
   init.Set("Month", prm.month);
   init.Set("Day", prm.day);
   init.Set("WeatherFilename", prm.fname);
   init.Set("Latitude", prm.lat);
   init.Set("Longitude", prm.lon);
   init.Set("AngstromIntercept", prm.angb0);
   init.Set("AngstromSlope", prm.angb1);

   // flux properties:
   init.Set("ReferenceHeight", prm.refhgt);
   init.Set("WindAttnCoefficient", prm.kwind);
   init.Set("EddyAttnCoefficient", prm.keddy);
   init.Set("StomataResA1", prm.rstA1);
   init.Set("StomataResA2", prm.rstA2);
   init.Set("MeanLeafWidth", prm.leafwidth);

   // soil properties:
   init.Set("PoreSizeDist", prm.poredist);
   init.Set("Porosity", prm.porosity);
   init.Set("WiltingPoint", prm.vwcwp);
   init.Set("SaturationPoint", prm.vwcsat);
   init.Set("WaterAmount01", prm.vwc01);
   init.Set("WaterAmount02", prm.vwc02);
   init.Set("Depth01", prm.len01);
   init.Set("Depth02", prm.len02);
   init.Set("MaxRootDepth", prm.maxlen);
   init.Set("HydraulicSlope", prm.alpha);
   init.Set("SaturatedHydraulic", prm.ksat);

   // plant properties:
   init.Set("DevStage", prm.dvs);
   init.Set("BaseTemp01", prm.basetemp01);
   init.Set("BaseTemp02", prm.basetemp02);
   init.Set("LAI", prm.lai);
   init.Set("Height", prm.hgt);
   init.Set("HeightMax", prm.hgtmax);
   init.Set("HeightIntercept", prm.hgtb0);
   init.Set("HeightSlope", prm.hgtb1);
   init.Set("RootsGrowRate", prm.rootgrate);
   init.Set("StemWgt", prm.dm.stem);
   init.Set("GreenLeavesWgt", prm.dm.gnleaves);
   init.Set("DeadLeavesWgt", prm.dm.ddleaves);
   init.Set("RootsWgt", prm.dm.roots);
   init.Set("StorageWgt", prm.dm.storage);

   // file names for tabulated data lookups:
   init.Set("DevRate01Table", prm.tabdvr01);
   init.Set("DevRate02Table", prm.tabdvr02);
   init.Set("StemTable", prm.tabst);
   init.Set("LeavesTable", prm.tablv);
   init.Set("RootsTable", prm.tabrt);
   init.Set("StorageTable", prm.tabso);
   init.Set("SLATable", prm.tabsla);

   return in;
}

#endif
