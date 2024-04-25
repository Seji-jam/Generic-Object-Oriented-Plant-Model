#ifndef INI_H
#define INI_H

#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "afgen.h"

// ini class
//    Reads variable names and set their values from input file
class ini
{
private:
   std::ifstream &in_;                 // data file
   std::ios_base::iostate oldstate_;   // data file original settings

   ini(ini const &);
   ini &operator=(ini const &);
   void Trim(std::string &str) const;
   std::string const Lowercase(std::string const &str) const;
   void Parse(std::string const &str, std::string &lhs,
              std::string &rhs, bool bString);
   void FindValue(std::string const &variable, std::string &value);

public:
   explicit ini(std::ifstream &file);
   ~ini();

   template <typename T> void Set(std::string const &variable,
                                  T &value);
   void Set(std::string const &variable, bool &value);
   void Set(std::string const &variable, std::string &value);
   void Set(std::string const &variable, afgen &value);
};

// Constructs object and set the I/O flags
inline ini::ini(std::ifstream &file)
   : in_(file), oldstate_(file.exceptions())
{
   in_.exceptions(std::ios_base::failbit | std::ios_base::badbit |
                  std::ios_base::eofbit);
}

inline ini::~ini()
{
   in_.exceptions(oldstate_);   // restore the I/O flags
}

// Reads the variable name and sets its value from file
template <typename T>
void ini::Set(std::string const &variable, T &value)
{
   std::string str;
   FindValue(variable, str);
   std::stringstream ss(str);
   if (!(ss >> value))
      throw std::runtime_error("Bad value read");
}

// For booleans: reads the variable name and sets its value from file
inline void ini::Set(std::string const &variable, bool &value)
{
   std::string str;
   FindValue(variable, str);
   std::stringstream ss(Lowercase(str));
   if (!(ss >> std::boolalpha >> value))
      throw std::runtime_error("Bad value read");
}

// For strings: reads the variable name and sets its value from file
inline void ini::Set(std::string const &variable, std::string &value)
{
   std::string str;
   while (str.empty())
   {
      char buf[255];
      in_.getline(buf, 255);
      str = buf;
      Trim(str);
   }

   std::string lhs, rhs;
   Parse(str, lhs, rhs, true);

   if (Lowercase(lhs) != Lowercase(variable))
      throw std::runtime_error("Variable name not found");

   while (rhs.empty())
   {
      char buf[255];
      in_.getline(buf, 255);
      rhs = buf;
      Trim(rhs);
   }

   value = rhs;
}

// For afgen class: reads the file name and sets the afgen object
inline void ini::Set(std::string const &variable, afgen &value)
{
   std::string fname;      // filename for tabulated data
   Set(variable, fname);   // get the file name now
   value.Filename(fname);  // set the afgen object
   value.MemoriseTable();  // read in all tabulated data
}

#endif
