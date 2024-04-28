#include "clck.h"

clck &clck::operator=(clck const &rhs)
{
   date_ = rhs.date_;
   time_ = rhs.time_;
   timezone_ = rhs.timezone_;
   dlstime_ = rhs.dlstime_;
   return *this;
}

bool clck::operator==(clck const &rhs) const
{
   return (date_.year==rhs.date_.year &&
           date_.month==rhs.date_.month &&
           date_.day==rhs.date_.day &&
           time_.hour==rhs.time_.hour &&
           time_.minute==rhs.time_.minute);
}

bool clck::operator>(clck const &rhs) const
{
   bool bEqYr = (date_.year==rhs.date_.year);
   bool bEqMt = (date_.month==rhs.date_.month);
   bool bEqDy = (date_.day==rhs.date_.day);
   bool bEqHr = (time_.hour==rhs.time_.hour);
   return (date_.year>rhs.date_.year ||
          (bEqYr && date_.month>rhs.date_.month) ||
          (bEqYr && bEqMt && date_.day>rhs.date_.day) ||
          (bEqYr && bEqMt && bEqDy &&
           time_.hour>rhs.time_.hour) ||
          (bEqYr && bEqMt && bEqDy && bEqHr &&
           time_.minute>rhs.time_.minute));
}

bool clck::operator<(clck const &rhs) const
{
   bool bEqYr = (date_.year==rhs.date_.year);
   bool bEqMt = (date_.month==rhs.date_.month);
   bool bEqDy = (date_.day==rhs.date_.day);
   bool bEqHr = (time_.hour==rhs.time_.hour);
   return (date_.year<rhs.date_.year ||
          (bEqYr && date_.month<rhs.date_.month) ||
          (bEqYr && bEqMt && date_.day<rhs.date_.day) ||
          (bEqYr && bEqMt && bEqDy &&
           time_.hour<rhs.time_.hour) ||
          (bEqYr && bEqMt && bEqDy && bEqHr &&
           time_.minute<rhs.time_.minute));
}

// INTERNAL USE: Check validity of given year (nYear),
//   month (nMonth) and day (nDay)
GDATE const clck::CheckDate(int nYear, int nMonth, int nDay) const
{
   if (nYear<1900 || nMonth<1 || nMonth>12 || nDay<1 ||
       nDay>MonthLen(nYear, nMonth))
      throw std::invalid_argument(
      "clck::CheckDate(int, int, int). Invalid date given.");

   return GDATE(nYear, nMonth, nDay);
}

// INTERNAL USE: Check validity of given hour (nHour)
//   and minute (nMinute)
GTIME const clck::CheckTime(int nHour, int nMinute) const
{
   if (nHour<0 || nHour>23 || nMinute<0 || nMinute>59)
      throw std::invalid_argument(
      "clck::CheckTime(int, int). Invalid time given.");

   return GTIME(nHour, nMinute);
}

// INTERNAL USE: The month length; that is, no. of days in a
//   given month (nMonth) and year (nYear)
int clck::DaysInMonth(int nYear, int nMonth) const
{
   if (nMonth<1 || nMonth>12)
      throw std::invalid_argument(
      "clck::DaysInMonth(int, int). Invalid month given.");

   return MonthLen(nYear, nMonth);
}

// INTERNAL USE: No. of days of a given date (dt) since Jan. 1
int clck::DayOfYear(GDATE const &dt) const
{
   int const CUMDAY[12] = {0, 31, 59, 90, 120, 151, 181, 212,
                           243, 273, 304, 334};
   int cDoy = CheckDate(dt.year, dt.month, dt.day).day
              + CUMDAY[dt.month-1];
   if (dt.month>2 && IsLeapYear(dt.year))
      ++cDoy;
   return cDoy;
}

// INTERNAL USE: Adjust the month so it lies within range
void clck::CorrectTheMonth(GDATE &dt) const
{
   if (dt.month > 12)
   {
      while (dt.month > 12)
      {
         dt.month -= 12;
         ++dt.year;
      }
   }
   else if (dt.month <= 0)
   {
      while (dt.month <= 0)
      {
         dt.month += 12;
         --dt.year;
      }
   }
}

// INTERNAL USE: Adjust the day so it lies within range
void clck::CorrectTheDay(GDATE &dt) const
{
   int monthlen = clck::DaysInMonth(dt.year, dt.month);
   if (dt.day > monthlen)
   {
      while (dt.day > monthlen)
      {
         dt.day -= monthlen;
         ++dt.month;
         CorrectTheMonth(dt);
         monthlen = clck::DaysInMonth(dt.year, dt.month);
      }
   }
   else if (dt.day <= 0)
   {
      while (dt.day <= 0)
      {
         --dt.month;
         CorrectTheMonth(dt);
         monthlen = clck::DaysInMonth(dt.year, dt.month);
         dt.day += monthlen;
      }
   }
}

// INTERNAL USE: Adjust the hour so it lies within range
void clck::CorrectTheHour(GDATE &dt, GTIME &tm) const
{
   if (tm.hour >= 24)
   {
      while (tm.hour >= 24)
      {
         tm.hour -= 24;
         ++dt.day;
         CorrectTheDay(dt);
      }
   }
   else if (tm.hour < 0)
   {
      while (tm.hour < 0)
      {
         tm.hour += 24;
         --dt.day;
         CorrectTheDay(dt);
      }
   }
}

// INTERNAL USE: Adjust the minute so it lies within range
void clck::CorrectTheMinute(GDATE &dt, GTIME &tm) const
{
   if (tm.minute >= 60)
   {
      while (tm.minute >= 60)
      {
         tm.minute -= 60;
         ++tm.hour;
         CorrectTheHour(dt, tm);
      }
   }
   else if (tm.minute < 0)
   {
      while (tm.minute < 0)
      {
         tm.minute += 60;
         --tm.hour;
         CorrectTheHour(dt, tm);
      }
   }
}

// Returns in days the difference between two clocks (current and ck)
double clck::Difference(clck const &ck) const
{
   bool bSmaller = (*this < ck); 
   // ensure current clck (*this) preceeds ck
   clck const &cc0 = (bSmaller ? *this : ck);
   clck const &cc1 = (bSmaller ? ck : *this);
   double edoy0 = cc0.ExactDayOfYear();
   double edoy1 = cc1.ExactDayOfYear();
   int year0 = cc0.Year(), year1 = cc1.Year();
   double diff = edoy1 - edoy0 + ((year0==year1) ?
                                  0 : cc0.DaysInYear());
   for (int year=year0+1; year<year1; ++year)
      diff += DaysInYear(year);
   return diff;
}

// Advance current clck by days (note: seconds are ignored).
//   A negative days argument will retreat the clck.
void clck::Advance(double days)
{
   GDATE dt(date_);  // work on a copy
   GTIME tm(time_);
   // convert days into mins and add it all into current min
   tm.minute += static_cast<int>(days * 1440.0);
   // partition total minutes into hour, day, month and year
   CorrectTheMinute(dt, tm);
   date_ = dt;
   time_ = tm;
}
