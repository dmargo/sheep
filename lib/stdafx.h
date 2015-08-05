#pragma once

#include <cassert>

/* OPTION: Deduplicate edges while loading the graph.
 * NOTE: This option does not work with distributed loading. */
//#define DDUP_GRAPH

/* OPTION: Save preorder weight for each vertex in the tree.
 * These weights are needed by some (non-default) partitioning models.
 * However, they consume sizeof(esize_t) bytes of memory per vertex. */
//define USE_PRE_WEIGHT


#define USE_LLAMA
//#define USE_SNAP

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

