#include <iostream>
#include <fstream>

std::string set_fname (char inputfile[]) {
  std::string fname;
  fname.append(inputfile);
  fname = fname.substr(0, fname.rfind("."));
  return fname;
}
