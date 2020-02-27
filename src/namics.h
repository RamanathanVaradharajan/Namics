
#ifndef NAMICSxH
#define NAMICSxH

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
//#include <f2c.h>
#include <iomanip>

using namespace std;
//these are our options.
typedef double Real;


//nvcc Open.cu -lm -lcuda -lcudart -llapack -lblas -lf2c -lcublas -arch=sm_20 -o box
//I.V. Ionova, E.A. Carter, "Error vector choice in direct inversion in the iterative subspace method, J. Compt. Chem. 17, 1836-1847, 1996.

#ifndef MAINxH //here define global variables.
extern Real* BlasResult;
extern string version;
extern Real e;
extern Real T;
extern Real k_B;
extern Real k_BT;
extern Real PIE;
extern int DEBUG_BREAK;
extern Real eps0;
extern bool debug;

//extern Real factor;
#endif

enum MoleculeType {monomer, linear, branched, dendrimer, comb, ring};
enum transfer {to_segment, to_cleng, to_teng, to_bm, reset};
enum EngineType {SCF, CLENG, MESODYN, TENG};

template<typename T>
  auto load_argument_value(vector<string> args, string argument, T t) -> decltype(t) {
    vector<string>::iterator position;
    position = find(args.begin(), args.end(), argument);
    if ( position != args.end() ) {
      ++position;
      if (position != args.end() && (*position)[0] != '-') {
        istringstream ss(*position);
        ss >> t;
        return t;
      }
    }
    else throw 1;
    return t;
  }


#endif
/*
#ifdef _WIN32
    //define lapack for Windows (32-bit and 64-bit, this part is common)
   #ifdef _WIN64
      //define lapack for Windows (64-bit only)
   #else
      //define lapack for Windows (32-bit only)
   #endif
#elif __APPLE__
    //define lapack for apple
    #if TARGET_IPHONE_SIMULATOR
         // iOS Simulator
    #elif TARGET_OS_IPHONE
        // iOS device
    #elif TARGET_OS_MAC
        // Other kinds of Mac OS
    #else
    #   error "Unknown Apple platform"
    #endif
#elif __linux
    //#include <lapacke.h>
    #include <clapack.h>
#elif __unix__ // all unices not caught above
    #include <clapack.h>

#elif defined(_POSIX_VERSION)
    // POSIX
#else
#   error "Unknown compiler"
#endif
*/
