#include "wthrfile.h"
#include <sstream>

wthrfile::wthrfile(std::string fname, clck const &now)
   : wthrfname_(fname), now_(now)
{}

wthrfile::wthrfile(wthrfile const &rhs)
   : daymet_(rhs.daymet_), wthrfname_(rhs.wthrfname_),
     now_(rhs.now_)
{}

wthrfile &wthrfile::operator=(wthrfile const &rhs)
{
   daymet_ = rhs.daymet_;
   wthrfname_ = rhs.wthrfname_;
   now_ = rhs.now_;
   return *this;
}

//  INTERNAL USE: find the 4-digit year in a weather file name
//    fname  : weather file name
//  Returns the starting index position of year in file name
std::string::size_type
wthrfile::FindYearPos(std::string const &fname) const
{
   std::string::size_type npos = fname.find_last_of("\\/");
   if (npos == std::string::npos)
      npos = 0;  // no slashes
   else
      ++npos;
   return npos;
}

//  INTERNAL USE: form either next or previous year's
//    weather file name
//    fname   : weather file name
//    yeardir : set to NEXT or PREV for next or previous
//              year's file, respectively
//    pYear   : year of weather file (to be set inside function)
void
wthrfile::FilenameChangeYear(std::string &fname,
                             DateDir yeardir, int *pYear) const
{
   std::string::size_type npos = FindYearPos(fname);
   std::stringstream token;
   token << fname.substr(npos, 4);
   int year = 0;
   token >> year;
   if (yeardir == NEXT)
      ++year;
   else
      --year;

   std::stringstream token2;
   token2 << year;
   fname.replace(npos, 4, token2.str());

   if (pYear)
      *pYear = year;
}

//  INTERNAL USE: interpolate between two given values
//    clk   : current date
//    clk0  : start date
//    clk1  : end date
//    val0  : start value
//    val1  : end value
double wthrfile::Interpolate(clck const &clk, clck const &clk0,
                             clck const &clk1,
                             double val0, double val1) const
{
   double doy = clk.ExactDayOfYear();
   double doy0 = clk0.ExactDayOfYear();
   double doy1 = clk1.ExactDayOfYear();
   int year0 = clk0.Year(), year1 = clk1.Year();
   double rangedoy = 0.0;
   if (year0 == year1)
      rangedoy = (doy1 > doy0) ? (doy1 - doy0) : (doy0 - doy1);
   else
      rangedoy = (year1 > year0) ? 
                 (doy1 + clk0.DaysInYear() - doy0) :
                 (doy0 + clk1.DaysInYear() - doy1);
   double diff = (doy > doy0) ? (doy - doy0) : (doy0 - doy);
   double v = (val1 > val0) ? (val1 - val0) : (val0 - val1);
   double u = v * diff / rangedoy;
   return ((val1 > val0) ? val0 + u : val0 - u);
}

//  INTERNAL USE: set daily weather data from file; called from
//    DailyMeteo()
//    dw : data structure for possible interpolation operation
//    returns TRUE if read was successful, else FALSE
bool wthrfile::SetDailyMeteo(METEO_INTPOL_<DAILY_METEO> &dw) const
{
   bool bFound = false;
   
   std::ifstream in(dw.fname.c_str());    // open weather file
   if (!in)
      return bFound; // file not found or failed to open;
                     // exit with failure

   std::ios_base::iostate oldEx = in.exceptions();
   in.exceptions(std::ios_base::badbit);

   // work with variables;
   //   easier than calling clck access functions repeatedly
   clck clk = dw.clk;
   int year1 = dw.clk1.Year();
   int month = clk.Month(), month1 = 0;
   int day = clk.Day(), day1 = 0, doy = clk.DayOfYear();

   // Look for weather data for current date (clk).
   //   There are 2 possibilities:
   //   1. Exact match   -- month and day in file matches exactly
   //                       the current date
   //   2. Interpolation -- no exact match with current date, so
   //                       interpolate between the two closest
   //                       matching date
   if (year1 < dw.clk0.Year())
   {
      // we are looking at past year's weather file for
      //   interpolation, so do not find exact match; just go to
      //   the last line in file because it is the closest match
      //   with current date.
      while (in)
         // go to last line in file
         in >> month1 >> day1
            >> dw.met1.avsun >> dw.met1.tmax >> dw.met1.tmin
            >> dw.met1.rh >> dw.met1.wind >> dw.met1.rain;

      // *second closest match (see below **)
      dw.clk1.Date(GDATE(year1, month1, day1));
      dw.met.avsun = Interpolate(clk, dw.clk0, dw.clk1,
                                 dw.met0.avsun, dw.met1.avsun);
      dw.met.tmax = Interpolate(clk, dw.clk0, dw.clk1,
                                dw.met0.tmax, dw.met1.tmax);
      dw.met.tmin = Interpolate(clk, dw.clk0, dw.clk1,
                                dw.met0.tmin, dw.met1.tmin);
      dw.met.rh = Interpolate(clk, dw.clk0, dw.clk1,
                              dw.met0.rh, dw.met1.rh);
      dw.met.wind = Interpolate(clk, dw.clk0, dw.clk1,
                                dw.met0.wind, dw.met1.wind);
      dw.met.rain = Interpolate(clk, dw.clk0, dw.clk1,
                                dw.met0.rain, dw.met1.rain);
      bFound = true;
      in.close();
   }
   else
   {
      // look for exact match, with possibility of interpolation
      bool bFirstline = true;    // first time reading file
      bool bBackward = false;    // no need to look at past year
                                 //   weather file (not yet)
      while (!bFound && in)
      {
         // current line in file
         in >> month1 >> day1
            >> dw.met1.avsun >> dw.met1.tmax >> dw.met1.tmin
            >> dw.met1.rh >> dw.met1.wind >> dw.met1.rain;

         dw.clk1.Date(GDATE(year1, month1, day1));
         int doy1 = dw.clk1.DayOfYear();

         if (bFirstline && dw.bOriginalFile && (month<month1 ||
             (month==month1 && day<day1)))
         {
            // date in first line of weather file is already past
            //   the current date; this means we need to look at
            //   past year file for interpolation
            bBackward = true;
            // **first closest match (see above *)
            dw.clk0 = dw.clk1;  
            dw.met0 = dw.met1;
            int yr = 0;
            FilenameChangeYear(dw.fname, PREV, &yr);
            dw.clk1.Date(GDATE(yr, 1, 1));
            break;
         }
         else
         {
            // exact match?
            if (dw.bOriginalFile && month==month1 && day==day1)
            {
               dw.met = dw.met1;
               bFound = true;
            }

            // if exact match not found, try to interpolate, but
            //   only do so:
            //   1) if at least two lines of the file has been
            //      read (!bFirstline) or
            //   2) if current file is the next year file
            //      (!bOriginalFile)
            if (!bFound && (!bFirstline || !dw.bOriginalFile))
            {
               int doy0 = dw.clk0.DayOfYear();
               if (!dw.bOriginalFile || (doy > doy0 &&
                   doy < doy1))
               {
                  dw.met.avsun = Interpolate(clk, dw.clk0,
                                             dw.clk1,
                                             dw.met0.avsun,
                                             dw.met1.avsun);
                  dw.met.tmax = Interpolate(clk, dw.clk0,
                                            dw.clk1,
                                            dw.met0.tmax,
                                            dw.met1.tmax);
                  dw.met.tmin = Interpolate(clk, dw.clk0,
                                            dw.clk1,
                                            dw.met0.tmin,
                                            dw.met1.tmin);
                  dw.met.rh = Interpolate(clk, dw.clk0, dw.clk1,
                                          dw.met0.rh,
                                          dw.met1.rh);
                  dw.met.wind = Interpolate(clk, dw.clk0,
                                            dw.clk1,
                                            dw.met0.wind,
                                            dw.met1.wind);
                  dw.met.rain = Interpolate(clk, dw.clk0,
                                            dw.clk1,
                                            dw.met0.rain,
                                            dw.met1.rain);
                  bFound = true;
               }
            }

            // no longer reading the first line of file
            bFirstline = false;
            dw.clk0.Date(GDATE(year1, month1, day1));
            dw.met0 = dw.met1;
         }
      }

      in.close();

      // if match still not found (!bFound), try next or
      //   previous year's weather file, but do this only once
      //   (bOriginalFile) -- avoid looking in more than 2 files
      if (!bFound && dw.bOriginalFile)
      {
         dw.bOriginalFile = false;
         if (!bBackward)
         {
            int yr = 0;
            FilenameChangeYear(dw.fname, NEXT, &yr);
            dw.clk1.Date(GDATE(yr, 1, 1));
         }
         bFound = SetDailyMeteo(dw);  // recursive call
      }
   }

   in.exceptions(oldEx);   // reset
   return bFound;
}

//  INTERNAL USE: change the daily weather file name
//    (wthrfname_) to match the current year
void wthrfile::UpdateDailyWthrFname()
{
   std::string::size_type npos = FindYearPos(wthrfname_);
   std::stringstream token;
   token << now_.Year();
   wthrfname_.replace(npos, 4, token.str());
}

//  Read the weather file for the daily meteorological properties
//    for the current date (NOTE: Set the clck attribute first
//    to set the current date before calling this function).
//  Returns:
//    TRUE if read was successful, else FALSE
//    Note: daymet_ attribute will be set after function return
bool wthrfile::SetDailyMeteo()
{
   // ensure filename format matches current date
   UpdateDailyWthrFname();
   clck now(now_);           // current date and time
   METEO_INTPOL_<DAILY_METEO> dw(wthrfname_, now, now, now);

   if (SetDailyMeteo(dw))
   {
      // no problems, so set the meteorological attribute
      daymet_ = dw.met;
      return true;
   }

   return false;  // weather file read error
}
