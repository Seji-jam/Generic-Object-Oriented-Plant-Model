#ifndef CLCK_H
#define CLCK_H

#include <stdexcept>

// GDATE
struct GDATE
{
   int year, month, day;
   GDATE(int yr=2000, int mth=1, int dy=1)
      : year(yr), month(mth), day(dy) {}
};

// GTIME
struct GTIME
{
   int hour, minute;
   GTIME(int hr=0, int mn=0) : hour(hr), minute(mn) {}
};

// clck
class clck
{
private:
   GDATE date_;
   GTIME time_;
   int timezone_;  // time zone hours relative to GMT/UTC (hour)
   int dlstime_;   // daylight savings time (hour)

   int MonthLen(int nYear, int nMonth) const
   {
      int const MONTHLEN[12] = {31, (IsLeapYear(nYear) ? 29 : 28),
                                31, 30, 31, 30, 31, 31, 30, 31,
                                30, 31};
      return MONTHLEN[nMonth-1];
   }

   bool IsLeapYear(int nYear) const
   {
      return ((nYear%4==0 && nYear%100!=0) || nYear%400==0);
   }
   
   double ExactDayOfYear(GDATE const &dt, GTIME const &tm) const
   {
      return (DayOfYear(dt)+(tm.hour/24.0)+(tm.minute/1440.0));
   }

   int DaysInYear(int nYear) const
   {
      return (IsLeapYear(nYear) ? 366 : 365);
   }

   int DaysInMonth(int nYear, int nMonth) const;
   int DayOfYear(GDATE const &dt) const;
   GDATE const CheckDate(int nYear, int nMonth, int nDay) const;
   GTIME const CheckTime(int nHour, int nMinute) const;
   void CorrectTheMonth(GDATE &dt) const;
   void CorrectTheDay(GDATE &dt) const;
   void CorrectTheHour(GDATE &dt, GTIME &tm) const;
   void CorrectTheMinute(GDATE &dt, GTIME &tm) const;

public:
   clck(int yr=2000, int mth=1, int day=1, int hr=0,
        int mn=0, int tz=0, int dls=0)
      : date_(CheckDate(yr, mth, day)), time_(CheckTime(hr, mn)),
        timezone_(tz), dlstime_(dls)
   {}

   clck(clck const &rhs)
      : date_(rhs.date_), time_(rhs.time_),
        timezone_(rhs.timezone_), dlstime_(rhs.dlstime_)
   {}

   clck &operator=(clck const &rhs);
   bool operator==(clck const &rhs) const;
   bool operator>(clck const &rhs) const;
   bool operator<(clck const &rhs) const;
   bool operator!=(clck const &rhs) const
   {
       return (!operator==(rhs));
   }
   bool operator>=(clck const &rhs) const
   {
      return (operator>(rhs) || operator==(rhs));
   }
   bool operator<=(clck const &rhs) const
   {
      return (operator<(rhs) || operator==(rhs));
   }

   GDATE const Date() const {return date_;}
   GTIME const Time() const {return time_;}
   int Year() const {return date_.year;}
   int Month() const {return date_.month;}
   int Day() const {return date_.day;}
   int Hour() const {return time_.hour;}
   int Minute() const {return time_.minute;}
   int TimeZone() const {return timezone_;}
   int DlsTime() const {return dlstime_;}

   void Date(GDATE const &dt)
   {
      date_ = CheckDate(dt.year, dt.month, dt.day);
   }
   void Time(GTIME const &tm)
   {
      time_ = CheckTime(tm.hour, tm.minute);
   }
   void TimeZone(int iTZ) {timezone_ = iTZ;}
   void DlsTime(int iDS) {dlstime_ = iDS;}

   bool IsLeapYear() const {return IsLeapYear(date_.year);}
   int DaysInYear() const  {return DaysInYear(date_.year);}
   int DaysInMonth() const
   {
      return DaysInMonth(date_.year, date_.month);
   }
   int DayOfYear() const {return DayOfYear(date_);}
   double ExactDayOfYear() const
   {
      return ExactDayOfYear(date_, time_);
   }
   double Difference(clck const &ck) const;
   void Advance(double days);
};

#endif
