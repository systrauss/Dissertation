#define timing_cxx
#include "timing.h"

//C++ libraries
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

//ROOT Libraries
#include "TString.h"

// class Timing
// {
//   public:
//     std::vector <std::vector < double > > geGe;
//     std::vector <std::vector < double > > geSiLi;
//     std::vector <std::vector < double > > siliSiLi;
//     Timing(int,int,char*);
// };

Timing::Timing ()
{
    //Just make a thing
}

Timing::Timing (int nGe, int nSiLi, char* cTiming)
{
    std::string buffer;
    int nCounter = 0;
    if(strcmp(cTiming,"")==0) //No file provided
    {
        cTiming = "Timing";
    }
    std::fstream fTiming(Form("user/%s.dat",cTiming)); //Timing File
    if(!fTiming.is_open())
    {
        std::cout << "Timing file did not open. Ending program." << std::endl;
        //Need to add a flag somehow that this is void.
    }
    getline (fTiming,buffer); //Label Line
    getline (fTiming,buffer); //Timing Line
    std::vector<double> row; //Row for adding a level in.
    while(strncmp(buffer.c_str(),"//ge-sili",2) != 0) //Loop through all the ge Detectors
    {
        geGe.push_back(row); //Put a new row of timing
        for(int k=0; k< 4; k++) //Loop through detector 1,  detector 2, low time, high time.
        {
            geGe[nCounter].push_back(std::atof(buffer.substr(0,buffer.find_first_of(',',0)).c_str())); //read in jth number in line
            buffer = buffer.substr(buffer.find_first_of(',',0)+1,buffer.find_first_of('\n',0)); //make a substring of the rest of the line
        }
        nCounter++;
        getline (fTiming,buffer); //Line with timing in it
    }
    std::cout << "Ge-ge timing read in" << std::endl;
    getline (fTiming,buffer); //Line with timing in it
    nCounter = 0;
    while(strncmp(buffer.c_str(),"//sili-sili",2) != 0)
    {
        geSiLi.push_back(row); //Put a new row of timing
        for(int k=0; k< 4; k++) //Loop through detector 1,  detector 2, low time, high time.
        {
            geSiLi[nCounter].push_back(std::atof(buffer.substr(0,buffer.find_first_of(',',0)).c_str())); //read in jth number in line
            buffer = buffer.substr(buffer.find_first_of(',',0)+1,buffer.find_first_of('\n',0)); //make a substring of the rest of the line
        }
        nCounter++;
        getline (fTiming,buffer); //Line with timing in it
    }
    std::cout << "Ge-sili timing read in" << std::endl;
    nCounter=0;
    getline (fTiming,buffer); //Line with timing in it
    while(!fTiming.eof()) //Loop through all the sili Detectors
    {
        siliSiLi.push_back(row); //Put a new row of timing
        for(int k=0; k< 4; k++) //Loop through detector 1,  detector 2, low time, high time.
        {
            siliSiLi[nCounter].push_back(std::atof(buffer.substr(0,buffer.find_first_of(',',0)).c_str())); //read in jth number in line
            buffer = buffer.substr(buffer.find_first_of(',',0)+1,buffer.find_first_of('\n',0)); //make a substring of the rest of the line
        }
        nCounter++;
        getline (fTiming,buffer); //Line with timing in it
    }
    std::cout << "Sili-sili timing read in" << std::endl;
}
