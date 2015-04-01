#pragma once

#include <algorithm>
#include <parallel/algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <limits>
#include <queue>
#include <unordered_set>
#include <vector>

#include "graph_wrapper.h"
#include "jnode.h"
#include "jtree.h"
#include "stdafx.h"

uint32_t simple_hash(uint32_t k) {
  return k % 2;
}

uint32_t knuth_hash(uint32_t k) {
  uint32_t prime = 2654435761;
  return k * prime;
}

uint32_t cormen_hash(uint32_t k) {
  double A = 0.5 * (sqrt(5) - 1);
  uint32_t s = floor(A * pow(2,32));
  return k * s;
}



class Partition {
public:
  short num_parts;
  std::vector<short> parts;
  #define INVALID_PART (short)-1

  inline Partition(std::vector<jnid_t> const &seq, JNodeTable &jnodes, short np,
      double balance_factor = 1.03, bool vtx_weight = false, bool pst_weight = true, bool pre_weight = false) :
    num_parts(np), parts(jnodes.size(), INVALID_PART)
  {
    size_t total_weight = 0;
    for (jnid_t id = 0; id != jnodes.size(); ++id) {
      if (vtx_weight) total_weight += 1;
      if (pst_weight) total_weight += jnodes.pst_weight(id);
      if (pre_weight) total_weight += jnodes.pre_weight(id);
    }
    size_t max_component = (total_weight / num_parts) * balance_factor;

    // For each jnid_t, assign a part.
    forwardPartition(jnodes, max_component, vtx_weight, pst_weight, pre_weight);
    //backwardPartition(jnodes, max_component, vtx_weight, pst_weight, pre_weight);

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



  // PARTITIONING ALGORITHMS
  inline void forwardPartition(JNodeTable &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    // Classic algorithm modified for FFD binpacking.
    // 1. Count the uncut component below X.
    // 2. If component_below(X) > max_component, pack bins.
    // Obviously there is some subtlety the bin packing setup. Things you might do:
    // 1. Pack kids more aggressively instead of halting as soon as the component fits (see batchtree.cpp)
    // 2. Try packing half-size components since this is ideal for bin packing.
    // 3. Move to edge-weighted stuff if you want to minimize edge cuts. This is probably the best idea.
    // 4. Spend more time reasoning about optimization criteria for communication volume.
    std::vector<size_t> part_size;
    std::vector<size_t> component_below(jnodes.size(), 0);
    for (jnid_t id = 0; id != jnodes.size(); ++id) {
      if (vtx_weight) component_below.at(id) += 1;
      if (pst_weight) component_below.at(id) += jnodes.pst_weight(id);
      if (pre_weight) component_below.at(id) += jnodes.pre_weight(id);

      if (component_below.at(id) > max_component) {
        std::sort(jnodes.kids(id).begin(), jnodes.kids(id).end(), [&](jnid_t const lhs, jnid_t const rhs)
          { return component_below.at(lhs) > component_below.at(rhs); });

        do {
          // Try to pack kids.
          for (auto itr = jnodes.kids(id).cbegin(); component_below.at(id) > max_component &&
                    itr!= jnodes.kids(id).cend(); ++itr)
          {
            jnid_t const kid = *itr;
            assert(component_below.at(kid) <= max_component);
            if (parts.at(kid) != INVALID_PART) continue;

            // Find a part (bin) for this kid.
            for (short cur_part = 0; cur_part != (short) part_size.size(); ++cur_part) {
              // If kid packs...
              if (part_size.at(cur_part) + component_below.at(kid) <= max_component) {
                component_below.at(id) -= component_below.at(kid);
                part_size.at(cur_part) += component_below.at(kid);
                parts.at(kid) = cur_part;
                break;
              }
            }
          }
          // If kid packing fails, open a new part (bin)
          if (component_below.at(id) > max_component)
            part_size.push_back(0);
        } while(component_below.at(id) > max_component);
      }
      assert(component_below.at(id) <= max_component);
      if (jnodes.parent(id) != INVALID_JNID)
        component_below.at(jnodes.parent(id)) += component_below.at(id);
    }

    for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
      if (parts.at(id) == INVALID_PART && jnodes.parent(id) != INVALID_JNID)
        parts.at(id) = parts.at(jnodes.parent(id));

      while (parts.at(id) == INVALID_PART) {
        for (short cur_part = part_size.size() - 1; cur_part != -1; --cur_part) {
          if (part_size.at(cur_part) + component_below.at(id) <= max_component) {
            part_size.at(cur_part) += component_below.at(id);
            parts.at(id) = cur_part;
            break;
          }
        }
        if (parts.at(id) == INVALID_PART)
          part_size.push_back(0);
      }
    }
  }

  inline void backwardPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    std::vector<size_t> component_below(jnodes.size(), 0);
    for (jnid_t id = 0; id != jnodes.size(); ++id) {
      if (vtx_weight) component_below.at(id) += 1;
      if (pst_weight) component_below.at(id) += jnodes.pst_weight(id);
      if (pre_weight) component_below.at(id) += jnodes.pre_weight(id);

      if (jnodes.parent(id) != INVALID_JNID)
        component_below.at(jnodes.parent(id)) += component_below.at(id);
    }

    jnid_t critical = std::distance(
      component_below.cbegin(),
      std::max_element(component_below.cbegin(), component_below.cend()));
    while (jnodes.kids(critical).size() != 0) {
      critical = *std::max_element(jnodes.kids(critical).cbegin(), jnodes.kids(critical).cend(),
        [&component_below](jnid_t const lhs, jnid_t const rhs) {
        return component_below.at(lhs) < component_below.at(rhs); });
      component_below.at(jnodes.parent(critical)) -= component_below.at(critical);
    }

    short cur_part = 0;
    size_t part_size = 0;
    while (critical != INVALID_JNID) {
      if (part_size + component_below.at(critical) < max_component) {
        parts.at(critical) = cur_part;
        part_size += component_below.at(critical);
      } else {
        parts.at(critical) = ++cur_part;
        part_size = component_below.at(critical);
      }
      critical = jnodes.parent(critical);
    }

    for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
      if (parts.at(id) == INVALID_PART)
        parts.at(id) = jnodes.parent(id) != INVALID_JNID ? parts.at(jnodes.parent(id)) : cur_part;
    }
  }

  //XXX This has been somewhat compelling for reducing CV; consider why.
  inline void depthPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    std::vector<size_t> depth(jnodes.size(), 0);
    for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id)
      if (jnodes.parent(id) != INVALID_JNID)
        depth.at(id) = depth.at(jnodes.parent(id)) + 1;

    std::vector<jnid_t> jnids(jnodes.size());
    std::iota(jnids.begin(), jnids.end(), 0);
    __gnu_parallel::stable_sort(jnids.begin(), jnids.end(),
      [&](jnid_t const lhs, jnid_t const rhs) { return depth.at(lhs) > depth.at(rhs); });

    short cur_part = 0;
    size_t cur_size = 0;
    for (size_t idx = 0; idx != jnids.size(); ++idx) {
      parts.at(jnids.at(idx)) = cur_part;
      if (vtx_weight) cur_size += 1;
      if (pst_weight) cur_size += jnodes.pst_weight(jnids.at(idx));
      if (pre_weight) cur_size += jnodes.pre_weight(jnids.at(idx));
      if (cur_size >= max_component) {
        ++cur_part;
        cur_size = 0;
      }
    }
  }

  //XXX This is practically anti-optimal...why?
  inline void heightPartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    std::vector<size_t> height(jnodes.size(), 0);
    for (jnid_t id = 0; id != jnodes.size(); ++id)
      if (jnodes.parent(id) != INVALID_JNID)
        height.at(jnodes.parent(id)) = std::max(height.at(jnodes.parent(id)), height.at(id) + 1);

    std::vector<jnid_t> jnids(jnodes.size());
    std::iota(jnids.begin(), jnids.end(), 0);
    __gnu_parallel::stable_sort(jnids.begin(), jnids.end(),
      [&](jnid_t const lhs, jnid_t const rhs) { return height.at(lhs) < height.at(rhs); });

    short cur_part = 0;
    size_t cur_size = 0;
    for (size_t idx = 0; idx != jnids.size(); ++idx) {
      parts.at(jnids.at(idx)) = cur_part;
      if (vtx_weight) cur_size += 1;
      if (pst_weight) cur_size += jnodes.pst_weight(jnids.at(idx));
      if (pre_weight) cur_size += jnodes.pre_weight(jnids.at(idx));
      if (++cur_size == max_component) {
        ++cur_part;
        cur_size = 0;
      }
    }
  }

  inline void naivePartition(JNodeTable const &jnodes, size_t const max_component,
      bool const vtx_weight, bool const pst_weight, bool const pre_weight)
  {
    short cur_part = 0;
    size_t cur_size = 0;
    for (jnid_t id = 0; id != jnodes.size(); ++id) {
      parts.at(id) = cur_part;
      if (vtx_weight) cur_size += 1;
      if (pst_weight) cur_size += jnodes.pst_weight(id);
      if (pre_weight) cur_size += jnodes.pre_weight(id);
      if (++cur_size == max_component) {
        ++cur_part;
        cur_size = 0;
      }
    }
  }

  inline void randomPartition(size_t vertex_count) {
    srand(time(nullptr));

    assert(parts.size() == 0);
    parts.resize(vertex_count);
    for (auto itr = parts.begin(); itr != parts.end(); ++itr)
      *itr = rand() % num_parts;
  }

  inline void readPartition(char const *filename) {
    std::ifstream stream(filename);
    short p;

    assert(parts.size() == 0);
    while (stream >> p)
      parts.push_back(p);
  }



  inline void print() const
  {
    short max_part = 0;
    size_t first_part = 0;
    size_t second_part = 0;

    for (short part : parts) {
      max_part = std::max(max_part, part);
      if (part == 0)
        first_part += 1;
      else if (part == 1)
        second_part += 1;
    }

    printf("Actually created %d partitions.\n", max_part + 1);
    printf("First two partition sizes: %zu and %zu\n", first_part, second_part);
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

    //size_t ECV_up = 0;
    size_t ECV_down = 0;
    std::vector<size_t> down_balance(*std::max_element(parts.cbegin(), parts.cend()) + 1, 0);

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

        //ECV_up_nbrs.insert((X_pos > Y_pos) ? X_part : Y_part);
        ECV_down_nbrs.insert((X_pos < Y_pos) ? X_part : Y_part);
        if (X_pos < Y_pos) down_balance.at(X_part) += 1;

      }
      //ECV_up += ECV_up_nbrs.size() - 1;
      ECV_down += ECV_down_nbrs.size() - 1;
    }

    size_t max_down_bal = *std::max_element(down_balance.cbegin(), down_balance.cend());

    //printf("ECV(up)  : %zu (%f%%)\n", ECV_up, (double) ECV_up / graph.getEdges());
    printf("ECV(down): %zu (%f%%)\n", ECV_down, (double) ECV_down / graph.getEdges());
    printf("  balance: %zu (%f%%)\n", max_down_bal, (double) max_down_bal / (graph.getEdges() / num_parts));
  }



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
    buf.weight = 0.0;
    std::ofstream stream(output_filename, std::ios::binary | std::ios::trunc);

    for (size_t i = 0; i < seq.size(); ++i) {
      vid_t const X = seq[i];
      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const Y = *eitr;
        size_t const j = pos.at(Y);

        if (i <= j) {
          buf.tail = i;
          buf.head = j;
          stream.write((char*)&buf, sizeof(xs1));
        }
      }
    }
  }

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



  template <typename GraphType>
  void fennel(GraphType const &graph, std::vector<vid_t> const &seq,
      size_t const max_component, bool const edge_balanced)
  {
    double const n = graph.getNodes();
    double const m = 2 * graph.getEdges(); // # of DIRECTED edges; getEdges() returns UNDIRECTED#.
    double const k = num_parts;

    double const y = 1.5;
    double const a = edge_balanced ?
      n * pow(k / m, y) : // From KDD14 paper.
      m * (pow(k, y - 1.0) / pow(n, y)); // From original FENNEL paper.

    std::vector<double> part_value;
    std::vector<double> part_size(num_parts, 0.0);

    //for (vid_t const X : seq) {
    for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
      vid_t const X = *nitr;
      double X_weight = edge_balanced ? ((double) graph.getDeg(X)) : 1.0;

      part_value.assign(num_parts, 0.0);
      for (auto X_itr = graph.getEdgeItr(X); !X_itr.isEnd(); ++X_itr) {
        vid_t const Y = *X_itr;
        //if (X <= Y) continue;
        if (parts[Y] != INVALID_PART)
          part_value[parts[Y]] += 1.0;
      }

      short max_part = 0;
      double max_value = std::numeric_limits<double>::lowest();
      for (short p = 0; p != num_parts; ++p) {
        if (part_size[p] + X_weight > max_component) continue; // Hard balance limit.

        //double p_cost = a * y * pow(part_size[p], y - 1.0); // From original FENNEL paper.
        double p_cost = a * pow(part_size[p] + X_weight, y) - a * pow(part_size[p], y);
        double p_value = part_value[p] - p_cost;
          if (p_value > max_value) {
          max_part = p;
          max_value = p_value;
        }

        if (part_size[p] == 0.0) break; // Everything will be 0.0 after this point.
      }
      parts[X] = max_part;
      part_size[max_part] += X_weight;
    }
  }

  void fennel(char const *const filename) {
    // I tried to privilege edge-partitioned fennel but hardcoding |V| and |E|
    // so that only one scan of the graph file would be necessary.
    // This is obviously cheesy, but it was a prototype, and even with this
    // privilege it proves too slow.
    vid_t max_vid = 4036529;
    vid_t vertex_count = 3997962;
    size_t edge_count = 34681189;
    double balance_factor = 1.03;

    size_t max_component = (edge_count/num_parts)*balance_factor;
    parts.assign(edge_count + 1, INVALID_PART);

    struct xs1 {
	    unsigned tail;
	    unsigned head;
	    float weight;
    };
    xs1 buf;
    /*
    std::ifstream stream(filename, std::ios::binary);
    while (!stream.eof()) {
      stream.read((char*)&buf, sizeof(xs1));
      if (buf.tail > max_vid)
        max_vid = buf.tail;
      if (buf.head > max_vid)
        max_vid = buf.head;
      ++edge_count;
    }
    */

    double const n = vertex_count;
    double const m = 2 * edge_count; // # of DIRECTED edges; getEdges() returns UNDIRECTED#.
    double const k = num_parts;

    double const y = 1.5;
    double const a = m * (pow(k, y - 1.0) / pow(n, y));

    std::vector<double> part_value;
    std::vector<double> part_size(num_parts, 0.0);
    std::vector<bool> touches_part(num_parts * (max_vid + 1));

    std::ifstream stream(filename, std::ios::binary);
    for (size_t eid = 0; !stream.eof(); ++eid) {
      stream.read((char*)&buf, sizeof(xs1));
      vid_t const X = buf.tail;
      vid_t const Y = buf.head;

      part_value.assign(num_parts, 0.0);
      for (short p = 0; k != num_parts; ++p) {
        if (touches_part.at(num_parts * X + p) == true)
          part_value[p] += 1.0;
        if (touches_part.at(num_parts * Y + p) == true)
          part_value[p] += 1.0;
      }

      short max_part = 0;
      double max_value = std::numeric_limits<double>::lowest();
      for (short p = 0; p != num_parts; ++p) {
        if (part_size[p] + 1.0 > max_component) continue; // Hard balance limit.

        double p_cost = a * pow(part_size[p] + 1.0, y) - a * pow(part_size[p], y);
        double p_value = part_value[p] - p_cost;
          if (p_value > max_value) {
          max_part = p;
          max_value = p_value;
        }

        if (part_size[p] == 0.0) break; // Everything will be 0.0 after this point.
      }

      parts[eid] = max_part;
      part_size[max_part] += 1.0;
      touches_part.at(num_parts * X + max_part) = true;
      touches_part.at(num_parts * X + max_part) = true;
    }
  }
};

