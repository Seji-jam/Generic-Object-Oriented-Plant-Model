#ifndef AFGEN_H
#define AFGEN_H

#include <string>
#include <map>

// look up tabulated data stored in a file or memory
typedef enum {TABLE_IN_FILE=0, TABLE_IN_MEMORY} TableSource;

class afgen
{
private:
   std::string fname_;     // file name where tabulated data is stored
   std::map<double, double> table_;  // stored tabulated data

   void Match(double const &key, double &key0, double &key1,
              double &val, double &val0, double &val1,
              bool &bFound, bool &bFirstLine) const;

public:
   afgen() {}
   explicit afgen(std::string const &name) : fname_(name) {}
   afgen(afgen const &rhs) : fname_(rhs.fname_), table_(rhs.table_) {}

   afgen &operator=(afgen const &rhs)
   {
      fname_ = rhs.fname_; table_ = rhs.table_; return *this;
   }

   std::string const &Filename() const {return fname_;}

   void Filename(std::string const &fname)
   {
      fname_ = fname; table_.clear();
   }

   std::map<double, double> const &Table() const {return table_;}

   void MemoriseTable();
   void ForgetTable() {table_.clear();}

   void LookUp(double const &key, double &val,
               TableSource source=TABLE_IN_MEMORY) const;

   bool operator()(double const &key, double &val,
                   TableSource source=TABLE_IN_MEMORY) const;
};

#endif
