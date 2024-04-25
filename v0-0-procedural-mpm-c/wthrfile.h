#ifndef WTHRFILE_H
#define WTHRFILE_H

#include <string>
#include <fstream>
#include "clck.h"

// DAILY_METEO
struct DAILY_METEO
{
   double avsun;        // sunshine hours (hour)
   double tmax, tmin;   // max. and min. air temperature (deg. C)
   double rh;           // relative humidity (%)
   double wind;         // wind speed (m/s)
   double rain;         // rainfall (mm)
   DAILY_METEO()
      : avsun(0.0), tmax(0.0), tmin(0.0), rh(0.0), wind(0.0) {}
};

// wthrfile
class wthrfile
{
private:
   // attributes:
   DAILY_METEO daymet_;     // daily meteorological properties
   std::string wthrfname_;  // daily weather file name
   clck now_;               // current date and time

   template <typename T>
   struct METEO_INTPOL_
   {
      std::string fname;
      clck clk, clk0, clk1;
      T met;
      T met0, met1;
      bool bOriginalFile;
      METEO_INTPOL_(std::string const &name, clck const &c,
                    clck const &c0, clck const &c1)
         : fname(name), clk(c), clk0(c0),
           clk1(c1), bOriginalFile(true) {}
   };

   double Interpolate(clck const &clk, clck const &clk0,
                      clck const &clk1,
                      double val0, double val1) const;
   bool SetDailyMeteo(METEO_INTPOL_<DAILY_METEO> &dw) const;
   std::string::size_type FindYearPos(std::string const &fname) const;
   typedef enum {NEXT=0, PREV} DateDir;
   void FilenameChangeYear(std::string &fname, DateDir yeardir,
                           int *pYear=0) const;
   void UpdateDailyWthrFname();

public:
   wthrfile() {}
   wthrfile(std::string fname, clck const &now);
   wthrfile(wthrfile const &rhs);
   wthrfile &operator=(wthrfile const &rhs);

   // daily weather properties:
   bool SetDailyMeteo();
   DAILY_METEO &DailyMeteo() {return daymet_;}
   DAILY_METEO const &DailyMeteo() const {return daymet_;}

   clck &Now() {return now_;}
   clck const &Now() const {return now_;}
   std::string const Filename() const {return wthrfname_;}
   void Filename(std::string fname) {wthrfname_ = fname;}
};

#endif
