#include <stdexcept>
#include <fstream>
#include "afgen.h"

// INTERNAL USE: Used by operator().
void afgen::Match(double const &key, double &key0, double &key1,
                  double &val, double &val0, double &val1,
                  bool &bFound, bool &bFirstLine) const
{
   val = val1;
   bFound = (key == key1);    // check if exact match
   // may interpolate if no exact match and
   //   the given key is between a smaller and larger key
   if (!bFound && !bFirstLine && (key > key0 && key < key1))
   {
      bFound = true;
      val = (key - key0) * (val1 - val0) / (key1 - key0);
      val += val0;
   }
   bFirstLine = false;
   key0 = key1;
   val0 = val1;
}

// Look up a value (val) based on a key (key) from tabulated data.
//   Will interpolate if exact key is not found. This function will
//   throw an exception if lookup/interpolation is not successful.
//   Use the overloaded operator() if throwing an exception is not
//   wanted. Set parameter source to TABLE_IN_FILE if the table is
//   to be read from file, or TABLE_IN_MEMORY if the table is to be
//   read from memory. If the latter, ensure MemoriseTable() has been
//   called prior to this function call.
void afgen::LookUp(double const &key, double &val,
                   TableSource source) const
{
   if (!operator()(key, val, source))
      throw std::runtime_error("Key not found in file.");
}

// Look up a value (val) based on a key (key) from tabulated data.
//   Will interpolate if exact key is not found. This function will
//   throw an exception if lookup/interpolation is not successful.
//   Use the function LookUp() if throwing an exception is required
//   for failed lookups. Set parameter source to TABLE_IN_FILE if
//   the table is to be read from file, or TABLE_IN_MEMORY if the
//   table is to be read from memory. If the latter, ensure
//   MemoriseTable() has been called prior to this function call.
bool afgen::operator()(double const &key, double &val,
                       TableSource source) const
{
   double key0 = double(), key1 = double();
   double val0 = double(), val1 = double();
   bool bFound = false, bFirstLine = true;

   if (source == TABLE_IN_MEMORY)
   {
      for (std::map<double, double>::const_iterator p=table_.begin();
           !bFound && p!=table_.end(); ++p)
      {
         key1 = p->first;
         val1 = p->second;
         Match(key, key0, key1, val, val0, val1, bFound, bFirstLine);
      }
   }
   else
   {
      std::ifstream in(fname_.c_str());
      if (!in)
         throw std::runtime_error("File open error.");

      while (!bFound && in)
      {
         in >> key1 >> val1;
         Match(key, key0, key1, val, val0, val1, bFound, bFirstLine);
      }

      in.close();
   }

   return bFound;
}

// Read tabulated data from file into map container for fast access.
void afgen::MemoriseTable()
{
   std::ifstream in(fname_.c_str());
   if (!in)
      throw std::runtime_error("File open error.");

   double key = double();
   double val = double();
   table_.clear();
   while (in)
   {
      in >> key >> val;
      table_.insert(std::make_pair(key, val));
   }
   in.close();
}
