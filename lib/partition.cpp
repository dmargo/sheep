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

#include "partition.h"

#include <cstdlib>
#include <limits>

#include <mpi.h>
#include <parallel/algorithm>

size_t get_weight(JNodeTable const &jnodes, jnid_t id,
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

Partition::Partition(std::vector<jnid_t> const &seq, JNodeTable &jnodes, part_t np,
    double balance_factor, bool vtx_weight, bool pst_weight, bool pre_weight) :
  parts(jnodes.size(), INVALID_PART), num_parts(np)
{
  size_t total_weight = 0;
  for (jnid_t id = 0; id != jnodes.size(); ++id)
    total_weight += get_weight(jnodes, id, vtx_weight, pst_weight, pre_weight);
  size_t max_component = (total_weight / num_parts) * balance_factor;

  // For each jnid_t, assign a part.
  forwardPartition(jnodes, max_component, vtx_weight, pst_weight, pre_weight);

  // Convert jnid_t-indexed parts to vid_t indexed parts.
  std::vector<part_t> tmp(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_PART);
  for (size_t i = 0; i != seq.size(); ++i)
    tmp.at(seq.at(i)) = parts.at(i);
  parts = std::move(tmp);
}

void Partition::mpi_sync()
{
  vid_t max_vid = parts.size();
  MPI_Datatype MPI_vid_t = sizeof(vid_t) == 4 ? MPI_UINT32_T : MPI_UINT64_T;
  MPI_Bcast(&max_vid, 1, MPI_vid_t, 0, MPI_COMM_WORLD);
  parts.resize(max_vid, INVALID_PART);

  assert(sizeof(part_t) == 2);
  MPI_Bcast(parts.data(), parts.size(), MPI_SHORT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&num_parts, 1, MPI_SHORT, 0, MPI_COMM_WORLD);
}



/* 
 * TREE PARTITIONING ALGORITHMS
 */
void Partition::forwardPartition(JNodeTable &jnodes, size_t const max_component,
    bool const vtx_weight, bool const pst_weight, bool const pre_weight)
{
  // Classic algorithm modified for FFD binpacking.
  // 1. Count the uncut component below X.
  // 2. If component_below(X) > max_component, pack bins.
  // Obviously there is some subtlety the bin packing setup. Things you might do:
  // 1. Pack kids more aggressively instead of halting as soon as the component fits.
  // 2. Try packing half-size components since this is ideal for bin packing.
  // 3. Move to edge-weighted stuff if you want to minimize edge cuts.
  // 4. Spend more time reasoning about optimization criteria for communication volume.
  std::vector<size_t> part_size;
  std::vector<size_t> component_below(jnodes.size(), 0);
  for (jnid_t id = 0; id != jnodes.size(); ++id)
  {
    component_below.at(id) += get_weight(jnodes, id, vtx_weight, pst_weight, pre_weight);
    if (component_below.at(id) > max_component)
    {
      std::sort(jnodes.kids(id).begin(), jnodes.kids(id).end(),
      [&component_below](jnid_t const lhs, jnid_t const rhs) {
        return component_below.at(lhs) > component_below.at(rhs); });

      // Try to pack kids
      do {
        for (auto itr = jnodes.kids(id).cbegin(); component_below.at(id) > max_component &&
                  itr!= jnodes.kids(id).cend(); ++itr)
        {
          jnid_t const kid = *itr;
          assert(component_below.at(kid) <= max_component);
          if (parts.at(kid) != INVALID_PART) continue;

          // Find a part (bin) for this kid.
          for (part_t cur_part = 0; cur_part != (part_t) part_size.size(); ++cur_part) {
            // If kid packs...
            if (part_size.at(cur_part) + component_below.at(kid) <= max_component) {
              component_below.at(id) -= component_below.at(kid);
              part_size.at(cur_part) += component_below.at(kid);
              parts.at(kid) = cur_part;
              break;
            }
          }
        }
        // If packing fails, open a new part (bin)
        if (component_below.at(id) > max_component)
          part_size.push_back(0);
      } while(component_below.at(id) > max_component);
    }
    assert(component_below.at(id) <= max_component);
    if (jnodes.parent(id) != INVALID_JNID)
      component_below.at(jnodes.parent(id)) += component_below.at(id);
  }

  /* At the conclusion of the loop, parts are only assigned to "cut" vertices.
   * So, push part assignments down the tree to vertices that don't yet have them. */
  for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
    if (parts.at(id) == INVALID_PART && jnodes.parent(id) != INVALID_JNID)
      parts.at(id) = parts.at(jnodes.parent(id));

    // If id is a root, then it needs to be packed into a bin.
    while (parts.at(id) == INVALID_PART) {
      for (part_t cur_part = part_size.size() - 1; cur_part != -1; --cur_part) {
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

void Partition::backwardPartition(JNodeTable const &jnodes, size_t const max_component,
    bool const vtx_weight, bool const pst_weight, bool const pre_weight)
{
  std::vector<size_t> component_below(jnodes.size(), 0);
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    component_below.at(id) += get_weight(jnodes, id, vtx_weight, pst_weight, pre_weight);
    if (jnodes.parent(id) != INVALID_JNID)
      component_below.at(jnodes.parent(id)) += component_below.at(id);
  }

  // Find the critical path
  jnid_t critical = std::distance(
    component_below.cbegin(),
    std::max_element(component_below.cbegin(), component_below.cend()));
  while (jnodes.kids(critical).size() != 0) {
    critical = *std::max_element(jnodes.kids(critical).cbegin(), jnodes.kids(critical).cend(),
    [&component_below](jnid_t const lhs, jnid_t const rhs) {
      return component_below.at(lhs) < component_below.at(rhs); });
    component_below.at(jnodes.parent(critical)) -= component_below.at(critical);
  }

  // Pack parts along the critical path
  part_t cur_part = 0;
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

  // Pack any unpacked parts, but this method is fundamentally broken for multiple components ATM.
  for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
    if (parts.at(id) == INVALID_PART)
      parts.at(id) = jnodes.parent(id) != INVALID_JNID ? parts.at(jnodes.parent(id)) : cur_part;
  }
}

//XXX This has been somewhat compelling for reducing CV on the cheap.
void Partition::depthPartition(JNodeTable const &jnodes, size_t const max_component,
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

  part_t cur_part = 0;
  size_t cur_size = 0;
  for (size_t idx = 0; idx != jnids.size(); ++idx) {
    parts.at(jnids.at(idx)) = cur_part;
    cur_size += get_weight(jnodes, jnids.at(idx), vtx_weight, pst_weight, pre_weight);
    if (cur_size >= max_component) {
      ++cur_part;
      cur_size = 0;
    }
  }
}

//XXX In contrast, this is practically anti-optimal, which is interesting.
void Partition::heightPartition(JNodeTable const &jnodes, size_t const max_component,
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

  part_t cur_part = 0;
  size_t cur_size = 0;
  for (size_t idx = 0; idx != jnids.size(); ++idx) {
    parts.at(jnids.at(idx)) = cur_part;
    cur_size += get_weight(jnodes, jnids.at(idx), vtx_weight, pst_weight, pre_weight);
    if (cur_size >= max_component) {
      ++cur_part;
      cur_size = 0;
    }
  }
}

void Partition::naivePartition(JNodeTable const &jnodes, size_t const max_component,
    bool const vtx_weight, bool const pst_weight, bool const pre_weight)
{
  part_t cur_part = 0;
  size_t cur_size = 0;
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    parts.at(id) = cur_part;
    cur_size += get_weight(jnodes, id, vtx_weight, pst_weight, pre_weight);
    if (cur_size >= max_component) {
      ++cur_part;
      cur_size = 0;
    }
  }
}

void Partition::randomPartition(size_t vertex_count) {
  srand(time(nullptr));

  assert(parts.size() == 0);
  parts.resize(vertex_count);
  for (auto itr = parts.begin(); itr != parts.end(); ++itr)
    *itr = rand() % num_parts;
}



/*
 * FENNEL IMPLEMENTATIONS
 */
template <typename GraphType>
void Partition::fennel(GraphType const &graph, std::vector<vid_t> const &seq,
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

    part_t max_part = 0;
    double max_value = std::numeric_limits<double>::lowest();
    for (part_t p = 0; p != num_parts; ++p) {
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

void Partition::fennel(char const *const filename) {
  // I tried to privilege edge-partitioned fennel by hardcoding |V| and |E|
  // so that only one scan of the graph file would be necessary.
  // This is obviously cheesy, but it was a prototype, and even with this
  // advantage it was too slow to be worthwhile for our evaluation.
  vid_t max_vid = 4036529;
  size_t vertex_count = 3997962;
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
    for (part_t p = 0; k != num_parts; ++p) {
      if (touches_part.at(num_parts * X + p) == true)
        part_value[p] += 1.0;
      if (touches_part.at(num_parts * Y + p) == true)
        part_value[p] += 1.0;
    }

    part_t max_part = 0;
    double max_value = std::numeric_limits<double>::lowest();
    for (part_t p = 0; p != num_parts; ++p) {
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



/*
 * EVALUATORS
 */
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
void Partition::evaluate(GraphType const &graph) const {
  size_t edges_cut = 0;
  size_t Vcom_vol = 0;
  size_t ECV_hash = 0;

  part_t max_part = *std::max_element(parts.cbegin(), parts.cend());
  std::vector<size_t> vertex_balance(max_part + 1, 0);
  std::vector<size_t> hash_balance(max_part + 1, 0);

  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
    vid_t const X = *nitr;
    part_t const X_part = parts.at(X);
    assert(X_part != INVALID_PART);
    vertex_balance.at(X_part) += 1;

    std::unordered_set<part_t> Vcom_vol_nbrs = {X_part};
    std::unordered_set<part_t> ECV_hash_nbrs = {};

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      part_t const Y_part = parts.at(Y);
      assert(Y_part != INVALID_PART);

      if (X < Y && X_part != Y_part) ++edges_cut;
      Vcom_vol_nbrs.insert(Y_part);

      part_t hash_part = cormen_hash(X) < cormen_hash(Y) ? X_part : Y_part;
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
void Partition::evaluate(GraphType const &graph, std::vector<vid_t> const &seq) const {
  evaluate(graph);

  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  size_t ECV_down = 0;
  size_t ECV_up = 0;

  part_t max_part = *std::max_element(parts.cbegin(), parts.cend()) + 1;
  std::vector<size_t> down_balance(max_part, 0);
  std::vector<size_t> up_balance(max_part, 0);

  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
    vid_t const X = *nitr;
    jnid_t const X_pos = pos.at(X);
    part_t const X_part = parts.at(X);
    assert(X_part != INVALID_PART);

    std::unordered_set<part_t> ECV_down_nbrs = {};
    std::unordered_set<part_t> ECV_up_nbrs = {};

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      jnid_t const Y_pos = pos.at(Y);
      part_t const Y_part = parts.at(Y);
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



/*
 * INPUT/OUTPUT
 */
template <typename GraphType, typename WriterType>
void Partition::writeIsomorphicGraph(
    GraphType const &graph, std::vector<vid_t> seq,
    char const *const output_filename) const
{
  std::stable_sort(seq.begin(), seq.end(), [&](vid_t const lhs, vid_t const rhs)
      { return parts.at(lhs) < parts.at(rhs); });

  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  WriterType writer(output_filename);
  for (jnid_t X_pos = 0; X_pos < seq.size(); ++X_pos) {
    vid_t const X = seq[X_pos];

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      jnid_t const Y_pos = pos.at(Y);

      // XXX GraphType is an undirected graph, so no need to write edges twice.
      if (X_pos < Y_pos)
        writer.write(X_pos,Y_pos);
    }
  }
}

template <typename ReaderType, typename WriterType>
void Partition::writeIsomorphicGraph_template(
    char const *const input_filename, std::vector<vid_t> seq,
    char const *const output_filename) const
{
  std::stable_sort(seq.begin(), seq.end(), [&](vid_t const lhs, vid_t const rhs)
      { return parts.at(lhs) < parts.at(rhs); });

  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  vid_t X,Y;
  ReaderType reader(input_filename);
  WriterType writer(output_filename);
  while(reader.read(X,Y)) {
    jnid_t X_pos = pos.at(X);
    jnid_t Y_pos = pos.at(Y);
    writer.write(X_pos,Y_pos);
  }
}

template <typename WriterType>
void Partition::writeIsomorphicGraph(
    char const *const input_filename, std::vector<vid_t> const &seq,
    char const *const output_filename) const
{
  if (strcmp(".dat", input_filename + strlen(input_filename) - 4) == 0)
    writeIsomorphicGraph_template<XS1Reader,WriterType>(input_filename, seq, output_filename);
  else
    writeIsomorphicGraph_template<SNAPReader,WriterType>(input_filename, seq, output_filename);
}

template <typename GraphType, typename WriterType>
void Partition::writePartitionedGraph(
    GraphType const &graph, std::vector<vid_t> const &seq,
    char const *const output_prefix) const
{
  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  part_t const max_part = *std::max_element(parts.cbegin(), parts.cend());
  assert(max_part < 10000);

  std::vector<WriterType*> writers;
  for (part_t p = 0; p != max_part + 1; ++p) {
    char *output_filename = (char*)malloc(strlen(output_prefix) + 5);
    sprintf(output_filename, "%s%04d", output_prefix, p);
    writers.emplace_back(new WriterType(output_filename));
    free(output_filename);
  }

  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
    vid_t const X = *nitr;
    jnid_t const X_pos = pos.at(X);
    part_t const X_part = parts.at(X);
    assert(X_part != INVALID_PART);

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      if (X >= Y) continue;
      //XXX GraphType is an undirected graph, so no need to write edges twice.

      jnid_t const Y_pos = pos.at(Y);
      part_t const Y_part = parts.at(Y);
      assert(Y_part != INVALID_PART);

      part_t edge_part = X_pos < Y_pos ? X_part : Y_part;
      writers.at(edge_part)->write(X,Y);
    }
  }

  for (WriterType *writer : writers)
    delete writer;
}

template <typename ReaderType, typename WriterType>
void Partition::writePartitionedGraph_template(
    char const *const input_filename, std::vector<vid_t> const &seq,
    char const *const output_prefix) const
{
  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  part_t const max_part = *std::max_element(parts.cbegin(), parts.cend());
  assert(max_part < 10000);

  std::vector<WriterType*> writers;
  for (part_t p = 0; p != max_part + 1; ++p) {
    char *output_filename = (char*)malloc(strlen(output_prefix) + 5);
    sprintf(output_filename, "%s%04d", output_prefix, p);
    writers.emplace_back(new WriterType(output_filename));
    free(output_filename);
  }

  vid_t X,Y;
  ReaderType reader(input_filename);
  while(reader.read(X,Y)) {
    jnid_t X_pos = pos.at(X);
    jnid_t Y_pos = pos.at(Y);

    part_t X_part = parts.at(X);
    part_t Y_part = parts.at(Y);

    assert(X_part != INVALID_PART);
    assert(Y_part != INVALID_PART);

    part_t edge_part = X_pos < Y_pos ? X_part : Y_part;
    writers.at(edge_part)->write(X,Y);
  }

  for (WriterType *writer : writers)
    delete writer;
}

template <typename WriterType>
void Partition::writePartitionedGraph(
    char const *const input_filename, std::vector<vid_t> const &seq,
    char const *const output_prefix) const
{
  if (strcmp(".dat", input_filename + strlen(input_filename) - 4) == 0)
    writePartitionedGraph_template<XS1Reader,WriterType>(input_filename, seq, output_prefix);
  else
    writePartitionedGraph_template<SNAPReader,WriterType>(input_filename, seq, output_prefix);
}

