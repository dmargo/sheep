#pragma once

#include <algorithm>
#include <parallel/algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <vector>

#include "graph_wrapper.h"
#include "jnode.h"
#include "stdafx.h"

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

  

  /* EVALUATION AND I/O */

  inline void print() const
  {
    short max_part = *std::max_element(parts.cbegin(), parts.cend()) + 1;
    size_t first_part = std::count(parts.cbegin(), parts.cend(), 0);
    size_t second_part = std::count(parts.cbegin(), parts.cend(), 1);

    printf("Actually created %d partitions.\n", max_part);
    printf("First two partition sizes: %zu and %zu\n", first_part, second_part);
  }

  static inline uint32_t simple_hash(vid_t k) {
    return k % 2;
  }

  static inline uint32_t knuth_hash(vid_t k) {
    uint32_t prime = 2654435761;
    return k * prime;
  }

  static inline uint32_t cormen_hash(vid_t k) {
    double A = 0.5 * (sqrt(5) - 1);
    uint32_t s = floor(A * pow(2,32));
    return k * s;
  }
  template <typename GraphType>
  inline void evaluate(GraphType const &graph) const {
    size_t edges_cut = 0;
    size_t Vcom_vol = 0;
    size_t ECV_hash = 0;

    short max_part = *std::max_element(parts.cbegin(), parts.cend());
    std::vector<size_t> vertex_balance(max_part + 1, 0);
    std::vector<size_t> hash_balance(max_part + 1, 0);

    for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
      vid_t const X = *nitr;
      short const X_part = parts.at(X);
      assert(X_part != INVALID_PART);
      vertex_balance.at(X_part) += 1;

      std::unordered_set<short> Vcom_vol_nbrs = {X_part};
      std::unordered_set<short> ECV_hash_nbrs = {};

      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const Y = *eitr;
        short const Y_part = parts.at(Y);
        assert(Y_part != INVALID_PART);

        if (X < Y && X_part != Y_part) ++edges_cut;
        Vcom_vol_nbrs.insert(Y_part);

        short hash_part = cormen_hash(X) < cormen_hash(Y) ? X_part : Y_part;
        ECV_hash_nbrs.insert(hash_part);
        if (X < Y) hash_balance.at(hash_part) += 1;

      }
      Vcom_vol += Vcom_vol_nbrs.size() - 1;
      ECV_hash += ECV_hash_nbrs.size() - 1;
    }

    size_t max_vertex_bal = *std::max_element(vertex_balance.cbegin(), vertex_balance.cend());
    size_t max_hash_bal = *std::max_element(hash_balance.cbegin(), hash_balance.cend());

    //XXX Remember graph.getEdges includes self-edges for some graphs.
    printf("edges cut: %zu (%f%%)\n", edges_cut, (double) edges_cut / graph.getEdges());
    printf("Vcom. vol: %zu (%f%%)\n", Vcom_vol, (double) Vcom_vol / graph.getEdges());
    printf("  balance: %zu (%f%%)\n", max_vertex_bal, (double) max_vertex_bal / (graph.getNodes() / num_parts));
    printf("ECV(hash): %zu (%f%%)\n", ECV_hash, (double) ECV_hash / graph.getEdges());
    printf("  balance: %zu (%f%%)\n", max_hash_bal, (double) max_hash_bal / (graph.getEdges() / num_parts));
  }

  template <typename GraphType>
  inline void evaluate(GraphType const &graph, std::vector<vid_t> const &seq) const {
    evaluate(graph);

    std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
    for (size_t i = 0; i != seq.size(); ++i)
      pos[seq[i]] = i;

    size_t ECV_down = 0;
    size_t ECV_up = 0;

    short max_part = *std::max_element(parts.cbegin(), parts.cend()) + 1;
    std::vector<size_t> down_balance(max_part, 0);
    std::vector<size_t> up_balance(max_part, 0);

    for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
      vid_t const X = *nitr;
      jnid_t const X_pos = pos.at(X);
      short const X_part = parts.at(X);
      assert(X_part != INVALID_PART);

      std::unordered_set<short> ECV_down_nbrs = {};
      std::unordered_set<short> ECV_up_nbrs = {};

      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const Y = *eitr;
        jnid_t const Y_pos = pos.at(Y);
        short const Y_part = parts.at(Y);
        assert(Y_part != INVALID_PART);

        ECV_down_nbrs.insert((X_pos < Y_pos) ? X_part : Y_part);
        ECV_up_nbrs.insert((X_pos > Y_pos) ? X_part : Y_part);
        if (X_pos < Y_pos) down_balance.at(X_part) += 1;
        if (X_pos > Y_pos) up_balance.at(X_part) += 1;
      }
      ECV_down += ECV_down_nbrs.size() - 1;
      ECV_up += ECV_up_nbrs.size() - 1;
    }

    size_t max_down_bal = *std::max_element(down_balance.cbegin(), down_balance.cend());
    size_t max_up_bal = *std::max_element(up_balance.cbegin(), up_balance.cend());

    printf("ECV(down): %zu (%f%%)\n", ECV_down, (double) ECV_down / graph.getEdges());
    printf("  balance: %zu (%f%%)\n", max_down_bal, (double) max_down_bal / (graph.getEdges() / num_parts));
    printf("ECV(up)  : %zu (%f%%)\n", ECV_up, (double) ECV_up / graph.getEdges());
    printf("  balance: %zu (%f%%)\n", max_up_bal, (double) max_up_bal / (graph.getEdges() / num_parts));

  }



  /* This write-out method reorders the graph such that if part[X] < part[Y] then X < Y.
   * It uses seq for tie-breaks.
   */
  template <typename GraphType>
  inline void writeIsomorphicGraph(
      GraphType const &graph, std::vector<vid_t> seq,
      char const *const output_filename) const
  {
    std::stable_sort(seq.begin(), seq.end(), [&](vid_t const lhs, vid_t const rhs)
        { return parts.at(lhs) < parts.at(rhs); });

    std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
    for (size_t i = 0; i != seq.size(); ++i)
      pos[seq[i]] = i;

    struct xs1 {
	    unsigned tail;
	    unsigned head;
	    float weight;
    };
    xs1 buf;
    buf.weight = 1.0;
    std::ofstream stream(output_filename, std::ios::binary | std::ios::trunc);

    for (size_t i = 0; i < seq.size(); ++i) {
      vid_t const X = seq[i];
      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const Y = *eitr;
        size_t const j = pos.at(Y);

        if (i >= j) {
          buf.tail = i;
          buf.head = j;
          stream.write((char*)&buf, sizeof(xs1));
        }
      }
    }
  }

  /* This write-out method simply writes each partition to a separate file.
   * It also isomorphs the graph according to seq, which is almost always desirable.
   */
  template <typename GraphType>
  inline void writePartitionedGraph(
      GraphType const &graph, std::vector<vid_t> const &seq,
      char const *const output_prefix) const
  {
    short const max_part = *std::max_element(parts.cbegin(), parts.cend());
    assert(max_part < 100);
    std::vector<std::ofstream*> output_streams;
    for (short p = 0; p != max_part + 1; ++p) {
      char *output_filename = (char*)malloc(strlen(output_prefix) + 3);
      sprintf(output_filename, "%s%02d", output_prefix, p);
      output_streams.emplace_back(new std::ofstream(output_filename, std::ios::trunc));
      free(output_filename);
    }
    
    std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
    for (size_t i = 0; i != seq.size(); ++i)
      pos[seq[i]] = i;

    for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
      vid_t const X = *nitr;
      jnid_t const X_pos = pos.at(X);
      short const X_part = parts.at(X);
      assert(X_part != INVALID_PART);

      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const Y = *eitr;
        if (X >= Y) continue;

        jnid_t const Y_pos = pos.at(Y);
        short const Y_part = parts.at(Y);
        assert(Y_part != INVALID_PART);

        short edge_part = X_pos < Y_pos ? X_part : Y_part;

        *output_streams.at(edge_part) << X << " " << Y << std::endl;
      }
    }

    for (std::ofstream *stream : output_streams)
      delete stream;
  }
};

