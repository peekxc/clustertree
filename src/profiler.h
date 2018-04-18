#ifndef PROFILER_H
#define PROFILER_H

#include <chrono>
#include <vector>
extern std::vector<double> benchmarks;
extern std::vector<std::string> benchmark_descrip;
extern std::chrono::time_point<std::chrono::high_resolution_clock> start;
extern std::chrono::time_point<std::chrono::high_resolution_clock> finish;

#define DO_PROFILING
#undef DO_PROFILING

#ifdef DO_PROFILING
  #define START_PROFILE() start = std::chrono::high_resolution_clock::now();
  #define END_PROFILE(x)                                         \
  finish = std::chrono::high_resolution_clock::now();            \
  Rprintf(x, std::chrono::duration<double, std::milli>(finish - start).count());

  #define ADD_BENCHMARK_DESCRIPTION(x, i) benchmark_descrip.at(i) = x;
  #define START_TIME() start = std::chrono::high_resolution_clock::now();
  #define END_TIME(x)                                            \
  finish = std::chrono::high_resolution_clock::now();            \
  benchmarks[x] = benchmarks[x] + std::chrono::duration<double, std::milli>(finish - start).count();

  #define START_NEW_TIME() auto new_start = std::chrono::high_resolution_clock::now();
  #define END_NEW_TIME(x)                                         \
  auto new_finish = std::chrono::high_resolution_clock::now();   \
  benchmarks[x] = benchmarks[x] + std::chrono::duration<double, std::milli>(new_finish - new_start).count();\

#else
  #define START_PROFILE()
  #define END_PROFILE(x)
  #define ADD_BENCHMARK_DESCRIPTION(x, i)
  #define START_TIME()
  #define END_TIME(x)
  #define START_NEW_TIME()
  #define END_NEW_TIME(x)
#endif

// Define whether to debug
#define DEBUG 1
#undef DEBUG

#ifdef DEBUG
#define R_OUT(x, ...) Rprintf(x, __VA_ARGS__);                 \
  FILE * pFile;                                                \
  pFile = fopen ("log.txt","a");                               \
  fprintf (pFile, x, __VA_ARGS__);                             \
  fclose (pFile);
#else
#define R_OUT(x, ...)
#endif


#endif