/*
 * Copyright (c) 2015
 *      The President and Fellows of Harvard College.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_set>
#include <vector>

#include "defs.h"
#include "graph_wrapper.h"
#include "jnode.h"
#include "readerwriter.h"

typedef short part_t;
#define INVALID_PART (part_t)-1

class Partition {
public:
  std::vector<part_t> parts;
  part_t num_parts;

  Partition(std::vector<jnid_t> const &seq, JNodeTable &jnodes, part_t np,
      double balance_factor = 1.03, bool vtx_weight = false, bool pst_weight = true, bool pre_weight = false);

  inline Partition() : parts(), num_parts() {}

  inline Partition(std::vector<jnid_t> const &seq, char const *filename) : parts(), num_parts()
  {
    readPartition(filename);
    num_parts = *std::max_element(parts.cbegin(), parts.cend());

    // Convert jnid_t-indexed parts to vid_t indexed parts.
    std::vector<part_t> tmp(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_PART);
    for (size_t i = 0; i != seq.size(); ++i)
      tmp.at(seq.at(i)) = parts.at(i);
    parts = std::move(tmp);
  }

  template <typename GraphType>
  inline Partition(GraphType const &graph, std::vector<vid_t> const &seq, part_t np, 
      double balance_factor = 1.03, bool edge_balanced = true) :
    parts(graph.getMaxVid() + 1, INVALID_PART), num_parts(np)
  {
    size_t total_weight = edge_balanced ? 2 * graph.getEdges() : graph.getNodes();
    size_t max_component = (total_weight / num_parts) * balance_factor;
    fennel(graph, seq, max_component, edge_balanced);
  }

  inline Partition(char const *filename, part_t np) :
    parts(), num_parts(np)
  {
    fennel(filename);
  }

  void mpi_sync();



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
    part_t p;

    assert(parts.size() == 0);
    while (stream >> p)
      parts.push_back(p);
  }

  template <typename GraphType>
  void fennel(GraphType const &graph, std::vector<vid_t> const &seq,
      size_t const max_component, bool const edge_balanced);

  void fennel(char const *const filename);



  /* EVALUATORS */

  inline void print() const
  {
    part_t max_part = *std::max_element(parts.cbegin(), parts.cend()) + 1;
    size_t first_part = std::count(parts.cbegin(), parts.cend(), 0);
    size_t second_part = std::count(parts.cbegin(), parts.cend(), 1);

    printf("Actually created %d partitions.\n", max_part);
    printf("First two partition sizes: %zu and %zu\n", first_part, second_part);
  }

  template <typename GraphType>
  void evaluate(GraphType const &graph) const; 

  template <typename GraphType>
  void evaluate(GraphType const &graph, std::vector<vid_t> const &seq) const;
  


  /* INPUT/OUTPUT 
   * This write-out method reorders the graph such that if part[X] < part[Y] then X < Y.
   * It uses seq for tie-breaks.
   */
  template <typename GraphType, typename WriterType = SNAPWriter>
  void writeIsomorphicGraph(
      GraphType const &graph, std::vector<vid_t> seq,
      char const *const output_filename) const;

  template <typename ReaderType, typename WriterType = SNAPWriter>
  void writeIsomorphicGraph_template(
      char const *const input_filename, std::vector<vid_t> seq,
      char const *const output_filename) const;

  template <typename WriterType = SNAPWriter>
  void writeIsomorphicGraph(
      char const *const input_filename, std::vector<vid_t> const &seq,
      char const *const output_filename) const;

  /* This write-out method simply writes each partition to a separate file.
   * It also isomorphs the graph according to seq, which is almost always desirable.
   */
  template <typename GraphType, typename WriterType = SNAPWriter>
  void writePartitionedGraph(
      GraphType const &graph, std::vector<vid_t> const &seq,
      char const *const output_prefix) const;

  template <typename ReaderType, typename WriterType = SNAPWriter>
  void writePartitionedGraph_template(
      char const *const input_filename, std::vector<vid_t> const &seq,
      char const *const output_prefix) const;

  template <typename WriterType = SNAPWriter>
  void writePartitionedGraph(
      char const *const input_filename, std::vector<vid_t> const &seq,
      char const *const output_prefix) const;
};

#include "partition.cpp"

