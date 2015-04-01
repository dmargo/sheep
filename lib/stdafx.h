#pragma once

#include <cassert>


#define USE_LLAMA
//#define USE_SNAP

//#define DDUP_GRAPH
#define DDUP_PST

#define USE_BUFIND
//#define USE_SUFIND

#define USE_PRE_WEIGHT


#ifdef USE_LLAMA
class LLAMAGraph;
typedef LLAMAGraph GraphWrapper;
typedef uint32_t vid_t;
typedef uint32_t esize_t;
#endif
#ifdef USE_SNAP
class SNAPGraph;
typedef SNAPGraph GraphWrapper;
typedef int vid_t;
typedef size_t esize_t;
#endif
#define INVALID_VID ((vid_t)-1)


#define KILO (1024)
#define MEGA (1024 * KILO)
#define GIGA (1024 * MEGA)

