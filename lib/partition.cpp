#include "partition.h"

#include <parallel/algorithm>
#include <cstdlib>
#include <limits>

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

    // Pack floating components
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
 * I/O AND EVALUATION
 */
template <typename GraphType>
void Partition::writeIsomorphicGraph(
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

      if (i < j) {
        buf.tail = i;
        buf.head = j;
        stream.write((char*)&buf, sizeof(xs1));
      }
    }
  }
}

template <typename GraphType>
void Partition::writePartitionedGraph(
    GraphType const &graph, std::vector<vid_t> const &seq,
    char const *const output_prefix) const
{
  part_t const max_part = *std::max_element(parts.cbegin(), parts.cend());
  assert(max_part < 100);
  std::vector<std::ofstream*> output_streams;
  for (part_t p = 0; p != max_part + 1; ++p) {
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
    part_t const X_part = parts.at(X);
    assert(X_part != INVALID_PART);

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      if (X >= Y) continue;

      jnid_t const Y_pos = pos.at(Y);
      part_t const Y_part = parts.at(Y);
      assert(Y_part != INVALID_PART);

      part_t edge_part = X_pos < Y_pos ? X_part : Y_part;

      *output_streams.at(edge_part) << X << " " << Y << std::endl;
    }
  }

  for (std::ofstream *stream : output_streams)
    delete stream;
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
  for (size_t i = 0; i != seq.size(); ++i)
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

