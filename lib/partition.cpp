#include "partition.h"

void Partition::forwardPartition(JNodeTable &jnodes, size_t const max_component,
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

    // Pack floating components
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

  short cur_part = 0;
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

  short cur_part = 0;
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
  short cur_part = 0;
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
template void Partition::fennel(GraphWrapper const &graph, std::vector<vid_t> const &seq,
    size_t const max_component, bool const edge_balanced);

void Partition::fennel(char const *const filename) {
  // I tried to privilege edge-partitioned fennel by hardcoding |V| and |E|
  // so that only one scan of the graph file would be necessary.
  // This is obviously cheesy, but it was a prototype, and even with this
  // advantage it was too slow to be worthwhile for our evaluation.
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

