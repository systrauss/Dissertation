//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan 11 12:37:54 2018 by ROOT version 5.34/19
// from TTree evtTree/133Ba_90deg_facingdown_trig(coin1(siliall+ge1+ge2) + pulser)
// found on file: rootFiles/Data/run1.root
//////////////////////////////////////////////////////////

#ifndef timing_h
#define timing_h

//C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class Timing
{
  public:
    std::vector <std::vector < double > > geGe;
    std::vector <std::vector < double > > geSiLi;
    std::vector <std::vector < double > > siliSiLi;
    Timing();
    Timing(int,int,char*);
};

#endif

#ifdef timing_cxx
#endif // #ifdef timing_cxx
