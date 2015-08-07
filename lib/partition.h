#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_set>
#include <vector>

#include "defs.h"
#include "graph_wrapper.h"
#include "jnode.h"

class Partition {
public:
  short num_parts;
  std::vector<short> parts;
  #define INVALID_PART (short)-1

  inline size_t get_weight(JNodeTable const &jnodes, jnid_t id,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    size_t result = 0;
    if (vtx_weight) result += 1;
    if (pst_weight) result += jnodes.pst_weight(id);
    if (pre_weight)
      for (jnid_t kid : jnodes.kids(id))
        result += jnodes.pre_weight(kid);
    return result;
  }

  inline Partition(std::vector<jnid_t> const &seq, JNodeTable &jnodes, short np,
      double balance_factor = 1.03, bool vtx_weight = false, bool pst_weight = true, bool pre_weight = false) :
    num_parts(np), parts(jnodes.size(), INVALID_PART)
  {
    size_t total_weight = 0;
    for (jnid_t id = 0; id != jnodes.size(); ++id)
      total_weight += get_weight(jnodes, id, vtx_weight, pst_weight, pre_weight);
    size_t max_component = (total_weight / num_parts) * balance_factor;

    // For each jnid_t, assign a part.
    forwardPartition(jnodes, max_component, vtx_weight, pst_weight, pre_weight);

    // Convert jnid_t-indexed parts to vid_t indexed parts.
    std::vector<short> tmp(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_PART);
    for (size_t i = 0; i != seq.size(); ++i)
      tmp.at(seq.at(i)) = parts.at(i);
    parts = std::move(tmp);
  }

  inline Partition(std::vector<jnid_t> const &seq, char const *filename) : num_parts(), parts()
  {
    readPartition(filename);
    num_parts = *std::max_element(parts.cbegin(), parts.cend());

    // Convert jnid_t-indexed parts to vid_t indexed parts.
    std::vector<short> tmp(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_PART);
    for (size_t i = 0; i != seq.size(); ++i)
      tmp.at(seq.at(i)) = parts.at(i);
    parts = std::move(tmp);
  }

  template <typename GraphType>
  inline Partition(GraphType const &graph, std::vector<vid_t> const &seq, short np, 
      double balance_factor = 1.03, bool edge_balanced = true) :
    num_parts(np), parts(graph.getMaxVid() + 1, INVALID_PART)
  {
    size_t total_weight = edge_balanced ? 2 * graph.getEdges() : graph.getNodes();
    size_t max_component = (total_weight / num_parts) * balance_factor;
    fennel(graph, seq, max_component, edge_balanced);
  }

  inline Partition(char const *filename, short np) :
    num_parts(np), parts()
  {
    fennel(filename);
  }



  /* PARTITIONING ALGORITHMS
   * The following algorithms are all intended for use with a JNodeTable.
   * forwardPartition is the best method and the method described in our paper;
   * the others are all experiments. */

  void forwardPartition(JNodeTable &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight);

  void backwardPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight);

  void depthPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight);

  void heightPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight);

  void naivePartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight);

  void randomPartition(size_t vertex_count);



  /* The remaining partitioning algorithms are NOT for JNodeTables.
   * readPartition does what it says on the can;
   * the other two are competing implementations of Fennel. */

  inline void readPartition(char const *filename) {
    std::ifstream stream(filename);
    short p;

    assert(parts.size() == 0);
    while (stream >> p)
      parts.push_back(p);
  }

  template <typename GraphType>
  void fennel(GraphType const &graph, std::vector<vid_t> const &seq,
      size_t const max_component, bool const edge_balanced);

  void fennel(char const *const filename);

  

  /* I/O AND EVALUATION */

  inline void print() const
  {
    short max_part = *std::max_element(parts.cbegin(), parts.cend()) + 1;
    size_t first_part = std::count(parts.cbegin(), parts.cend(), 0);
    size_t second_part = std::count(parts.cbegin(), parts.cend(), 1);

    printf("Actually created %d partitions.\n", max_part);
    printf("First two partition sizes: %zu and %zu\n", first_part, second_part);
  }

  /* This write-out method reorders the graph such that if part[X] < part[Y] then X < Y.
   * It uses seq for tie-breaks.
   */
  template <typename GraphType>
  inline void writeIsomorphicGraph(
      GraphType const &graph, std::vector<vid_t> seq,
      char const *const output_filename) const;

  /* This write-out method simply writes each partition to a separate file.
   * It also isomorphs the graph according to seq, which is almost always desirable.
   */
  template <typename GraphType>
  inline void writePartitionedGraph(
      GraphType const &graph, std::vector<vid_t> const &seq,
      char const *const output_prefix) const;

  template <typename GraphType>
  inline void evaluate(GraphType const &graph) const; 

  template <typename GraphType>
  inline void evaluate(GraphType const &graph, std::vector<vid_t> const &seq) const;
};

#include "partition.cpp"

