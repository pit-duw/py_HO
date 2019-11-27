#include <iostream>
#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t")
{
  const std::size_t strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return ""; // no content

  const std::size_t strEnd = str.find_last_not_of(whitespace);
  const std::size_t strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

std::string reduce(const std::string& str,
                   const std::string& fill = " ",
                   const std::string& whitespace = " \t")
{
  // trim first
  std::string result = trim(str, whitespace);

  // replace sub ranges
  std::size_t beginSpace = result.find_first_of(whitespace);
  while (beginSpace != std::string::npos)
    {
      const std::size_t endSpace = result.find_first_not_of(whitespace, beginSpace);
      const std::size_t range = endSpace - beginSpace;

      result.replace(beginSpace, range, fill);

      const std::size_t newStart = beginSpace + fill.length();
      beginSpace = result.find_first_of(whitespace, newStart);
    }

  return result;
}
