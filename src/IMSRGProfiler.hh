#ifndef IMSRGProfiler_h
#define IMSRGProfiler_h 1

#include <unistd.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <map>

/// Profiling class with a global instance called profiler.
/// This is for keeping track of timing and memory usage, etc.
/// IMSRGProfiler methods should probably not be called inside parallel blocks.

using namespace std;

class IMSRGProfiler
{
 public:
  // timer and counter are declaired as static so that there's only one copy of each of them
  static map<string, double> timer; ///< For keeping timing information for various method calls
  static map<string, int> counter;

  map<string,size_t> CheckMem();
  void PrintTimes();
  void PrintCounters();
  void PrintAll();
  size_t MaxMemUsage();
};

#endif
