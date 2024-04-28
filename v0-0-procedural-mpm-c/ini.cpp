#include <cstring>
#include <cstdlib>
#include <cctype>
#include "ini.h"

// INTERNAL USE: Removes whitespaces, quotes and equal sign from
//   given string (str)
void ini::Trim(std::string &str) const
{
   while (str[0] == ' ' || str[0] == '=' || str[0] == '\t'
          || str[0] == '"')
      str.erase(0, 1);

   int cLen = static_cast<int>(str.length());
   while (cLen > 0 && (str[cLen-1] == ' ' || str[cLen-1] == '\t'
                       || str[cLen-1] == '"'))
      str.erase(--cLen, 1);
}

// INTERNAL USE: Converts and returns a string, all chars in lowercase
std::string const ini::Lowercase(std::string const &str) const
{
   std::string locase(str);
   for (std::string::iterator pos=locase.begin();
        pos!=locase.end(); ++pos)
      *pos = static_cast<char>(std::tolower(*pos));
   return locase;
};

// INTERNAL USE: Splits a given string (str) into two: left- (lhs)
//   and right hand side (rhs). Where to split the two is where
//   the equal sign (=) is located in str.
//   If bString is set to TRUE, the rhs string can contain spaces.
void ini::Parse(std::string const &str, std::string &lhs,
                std::string &rhs, bool bString)
{
   lhs.clear();      // ensure all empty
   rhs.clear();

   char buf[255];
   std::strcpy(buf, str.c_str());

   char *pToken = std::strtok(buf, " =\t");
   int nCnt = 0;
   while (pToken && nCnt<2)
   {
      if (++nCnt == 1)
         lhs = pToken;
      else
      {
         rhs = pToken;
         Trim(rhs);
         if (rhs.empty())
            --nCnt;
      }

      if (!bString)
         pToken = std::strtok(0, " =\t");
      else
         pToken = std::strtok(0, "\n");
   }
}

// INTERNAL USE: Returns the variable's value as string
void ini::FindValue(std::string const &variable, std::string &value)
{
   std::string str;
   in_ >> str;

   std::string lhs, rhs;
   Parse(str, lhs, rhs, false);

   if (Lowercase(lhs) != Lowercase(variable))
      throw std::runtime_error("Variable name not found");

   while (rhs.empty())
   {
      in_ >> lhs;
      std::string unused;
      Parse(lhs, rhs, unused, false);
   }

   value = rhs;
}

