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

#include "defs.h"

#ifdef USE_LLAMA
#include <llama.h>
class LLAMAGraph {
private:
  ll_database database;
  ll_mlcsr_ro_graph *G; //XXX This should be const but LLAMA's methods need const markup
  size_t num_nodes;

public:
  LLAMAGraph(char const *filename, size_t const part = 0, size_t const num_parts = 0,
             bool const is_undirected = true) : database("./db/"), num_nodes(0)
  {
    ll_loader_config loader_config;
    loader_config.lc_no_properties = true;
    loader_config.lc_partial_load_part = part;
    loader_config.lc_partial_load_num_parts = num_parts;
    if (is_undirected)
      loader_config.lc_direction = LL_L_UNDIRECTED_DOUBLE;
    #ifdef DDUP_GRAPH
      loader_config.lc_deduplicate = true;
    #endif

    ll_file_loaders loaders;
    ll_file_loader* loader = loaders.loader_for(filename);
    loader->load_direct(database.graph(), filename, &loader_config);

    G = &database.graph()->ro_graph();
    for (node_t n = 0; n != G->max_nodes(); ++n)
      if (G->out_degree(n) != 0) ++num_nodes;
  }

  LLAMAGraph(LLAMAGraph &&other) = delete;
  LLAMAGraph(LLAMAGraph const &other) = delete;

  LLAMAGraph& operator=(LLAMAGraph &&other) = delete;
  LLAMAGraph& operator=(LLAMAGraph const &other) = delete;

  inline vid_t getMaxVid() const {
    return G->max_nodes();
  }

  inline size_t getNodes() const {
    return num_nodes;
  }

  inline size_t getEdges() const {
    return G->max_edges(0) / 2;
  }

  inline bool isNode(vid_t X) const {
    return (node_t) X < G->max_nodes() && G->out_degree(X) != 0;
  }

  inline size_t getDeg(vid_t X) const {
    return G->out_degree(X);   
  }

  class NodeItr {
  private:
    ll_mlcsr_ro_graph *const G;
    node_t n;

  public:
    NodeItr(ll_mlcsr_ro_graph *graph) : G(graph), n(0) {
      while (G->out_degree(n) == 0 && n != G->max_nodes())
        n++;
    }

    inline vid_t operator*() const {
      return n;
    }

    inline vid_t operator++() {
      do {
        ++n;
      } while (G->out_degree(n) == 0 && n != G->max_nodes());
      return operator*();
    }

    inline vid_t operator++(int) {
      vid_t result = operator*();
      operator++();
      return result;
    }

    inline bool isEnd() const {
      return n == G->max_nodes();
    }
  };

  inline NodeItr getNodeItr() const {
    return NodeItr(G);
  }

  class EdgeItr {
  private:
    ll_mlcsr_ro_graph *const G;
    ll_edge_iterator edge_itr;
    edge_t e;

  public:
    EdgeItr(ll_mlcsr_ro_graph *graph, vid_t X) : G(graph) {
      G->out_iter_begin(edge_itr, X);
      e = G->out_iter_next(edge_itr);
    }

    inline vid_t operator*() const {
      return edge_itr.last_node;
    }

    inline vid_t operator++() {
      e = G->out_iter_next(edge_itr);
      return operator*();
    }

    inline vid_t operator++(int) {
      vid_t result = operator*();
      operator++();
      return result;
    }

    inline bool isEnd() const {
      return e == LL_NIL_EDGE;
    }
  };

  inline EdgeItr getEdgeItr(vid_t X) const {
    return EdgeItr(G, X);
  }
};
typedef LLAMAGraph GraphWrapper;


#elif USE_SNAP
#include <Snap.h>
#undef max
#undef min
class SNAPGraph {
private:
  PUNGraph G;

public:
  SNAPGraph(char const *filename) {
    G = TSnap::LoadEdgeList<PUNGraph>(filename);
  }

  SNAPGraph(SNAPGraph &&other) = delete;
  SNAPGraph(SNAPGraph const &other) = delete;

  SNAPGraph& operator=(SNAPGraph &&other) = delete;
  SNAPGraph& operator=(SNAPGraph const &other) = delete;

  inline vid_t getMaxVid() const {
    vid_t max_vid = 0;
    for (auto nitr = getNodeItr(); !nitr.isEnd(); ++nitr)
      max_vid = std::max(max_vid, *nitr);
    return max_vid;
  }

  inline size_t getNodes() const {
    return (size_t)G->GetNodes();
  }

  inline size_t getEdges() const {
    return (size_t)G->GetEdges();
  }

  inline bool isNode(vid_t X) const {
    return G->IsNode(X);
  }

  inline size_t getDeg(vid_t X) const {
    return G->GetNI(X).GetDeg();
  }

  class NodeItr {
  private:
    TUNGraph::TNodeI node_itr;
    TUNGraph::TNodeI const end_itr;

  public:
    NodeItr(PUNGraph G) : node_itr(G->BegNI()), end_itr(G->EndNI()) {}

    inline vid_t operator*() const {
      return node_itr.GetId();
    }

    inline vid_t operator++() {
      node_itr++;
      return operator*();
    }

    inline vid_t operator++(int) {
      vid_t result = operator*();
      operator++();
      return result;
    }

    inline bool isEnd() const {
      return node_itr == end_itr;
    }
  };

  inline NodeItr getNodeItr() const {
    return NodeItr(G);
  }

  class EdgeItr {
  private:
    TUNGraph::TNodeI const edge_itr;
    int i;

  public:
    EdgeItr(PUNGraph G, vid_t X) : edge_itr(G->GetNI(X)), i(0) {}

    inline vid_t operator*() const {
      return edge_itr.GetOutNId(i);
    }

    inline vid_t& operator++() {
      i++;
      return operator*();
    }

    inline vid_t operator++(int) {
      vid_t result = operator*();
      operator++();
      return result;
    }

    inline bool isEnd() const {
      return i == edge_itr.GetDeg();
    }
  };

  inline EdgeItr getEdgeItr(vid_t X) const {
    return EdgeItr(G, X);
  }
};
typedef SNAPGraph GraphWrapper;
#endif

