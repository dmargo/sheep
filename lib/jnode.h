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

#include <fstream>
#include <vector>

#include <mpi.h>

#include "defs.h"
#include "jdata.h"
#include "merge.h"
#include "unionfind.h"

typedef vid_t jnid_t;
#define INVALID_JNID ((jnid_t)-1)

#ifdef USE_SIMPLE_UF
typedef SimpleUnionFind<jnid_t> UnionFind;
#else
typedef FastUnionFind<jnid_t> UnionFind;
#endif

class JNodeTable {
private:
  enum class State { ALLOCATED, MAPPED, TEMPORARY };
  State nodes_state;

  struct JNode {
    jnid_t parent;
    esize_t pst_weight;
    #ifdef USE_PRE_WEIGHT 
    esize_t pre_weight;
    #endif
    
    inline JNode() :
      parent(INVALID_JNID), pst_weight(0) {
        #ifdef USE_PRE_WEIGHT
          pre_weight = 0;
        #endif
      }
  }; 
  jnid_t end_id;
  jnid_t max_id;
  JNode *nodes;

  JDataTable<jnid_t> kid_data;
  JDataTable<vid_t>  pst_data;
  JDataTable<vid_t>  jxn_data;
  
  //XXX Roots are non-essential so consider a lazy construction -- or better, move this into the algorithm.
  UnionFind roots;

public:
  /* CONSTRUCTORS */
  JNodeTable() = delete;
  JNodeTable(jnid_t max_jnids, bool init_kids, size_t memory_limit);
  JNodeTable(char const *filename, jnid_t max_jnids, bool init_kids, size_t memory_limit);
  JNodeTable(char const *filename);
  JNodeTable(JNode *n, jnid_t end);

  JNodeTable(JNodeTable &&other);
  JNodeTable(JNodeTable const &other);
  JNodeTable(JNodeTable const &other, jnid_t partial_end);
  JNodeTable& operator=(JNodeTable &&other) = delete;
  JNodeTable& operator=(JNodeTable const &other) = delete;
  
  ~JNodeTable();


  void save(char const *filename);
  void merge(JNodeTable const &lhs, JNodeTable const &rhs, bool make_kids = false);
  void mpi_merge(bool make_kids = false);
  template <bool make_kids>
  friend void mpi_merge_reduction(void *in, void *inout, int *len, MPI_Datatype *datatype);


  inline jnid_t size() const { return end_id; }

  inline jnid_t newJNode() {
    assert(0 <= end_id && end_id < max_id);
    new(&nodes[end_id]) JNode();
    return end_id++;
  }

  //XXX UnionFind can't be revoked, so a JNode should never be deleted if adopt() has been called.
  inline void deleteJNode(jnid_t id) {
    assert(id == end_id - 1);

    if (id < kid_data.size()) {
      for (jnid_t const kid : kids(id)) {
        parent(kid) = INVALID_JNID;
        pre_weight(kid) = 0;
      }
      kid_data.deleteJData(id);}
    if (id < pst_data.size())
      pst_data.deleteJData(id);
    if (id < jxn_data.size())
      jxn_data.deleteJData(id);

    end_id--;
  }
 

  /* JNODE SET/GET/CONSTRUCTION HELPERS */
  inline jnid_t& parent(jnid_t id) { return nodes[id].parent; }
  inline jnid_t const & parent(jnid_t id) const { return nodes[id].parent; }

  inline esize_t& pst_weight(jnid_t id) { return nodes[id].pst_weight; }
  inline esize_t const & pst_weight(jnid_t id) const { return nodes[id].pst_weight; }

  inline esize_t& pre_weight(jnid_t id) {
    #ifdef USE_PRE_WEIGHT
      return nodes[id].pre_weight;
    #else
      esize_t static faux_weight = 0;
      return faux_weight;
    #endif
  }

  inline esize_t const & pre_weight(jnid_t id) const {
    #ifdef USE_PRE_WEIGHT
      return nodes[id].pre_weight;
    #else
      esize_t static const faux_weight = 0;
      return faux_weight;
    #endif
  }

  //XXX UnionFind can't be revoked, so a JNode should never be deleted after adopt()
  inline void adopt(jnid_t kid, jnid_t const id) {
    kid = roots.unify(kid, id);
    if (kid != id)
      parent(kid) = id;
  }
  

  /* JDATA TABLE WRAPPERS */
  inline void newKids(jnid_t id, size_t max_size) { 
    size_t tmp = kid_data.newJData(max_size, false);
    assert(tmp == id);
  }

  inline JData<jnid_t>& kids(jnid_t id) { return kid_data[id]; }
  inline JData<jnid_t> const & kids(jnid_t id) const { return kid_data[id]; }

  inline void meetKid(jnid_t kid, jnid_t const id, size_t weight) {
    kid = roots.find(kid);
    pre_weight(kid) += weight;

    if (parent(kid) != id) {
      parent(kid) = id;
      kids(id).push_back(kid);
    }
  }

  inline void adoptKids(jnid_t const id) {
    kid_data.shrinkJData(id);
    for (jnid_t const kid : kids(id))
      roots.unify(kid, id);
  }

  inline void makeKids() {
    kid_data = std::move(JDataTable<jnid_t>(max_id));

    std::vector<size_t> kids_size(size(), 0);
    for (jnid_t id = 0; id != size(); ++id) {
      newKids(id, kids_size[id]);
      if (parent(id) != INVALID_JNID)
        kids_size.at(parent(id)) += 1;
    }

    for (jnid_t id = 0; id != size(); ++id) {
      if (parent(id) != INVALID_JNID)
        kids(parent(id)).push_back(id);
    }
  }


  inline void newPst(jnid_t id, size_t max_size) {
    size_t tmp = pst_data.newJData(max_size, true);
    assert(tmp == id);
  }
  
  inline JData<vid_t>& pst(jnid_t id) { return pst_data[id]; }
  inline JData<vid_t> const & pst(jnid_t id) const { return pst_data[id]; }

  inline void cleanPst(jnid_t id) {
    std::sort(pst(id).begin(), pst(id).end());
    #ifndef DDUP_GRAPH
      auto end = std::unique(pst(id).begin(), pst(id).end());
      pst(id).len = std::distance(pst(id).begin(), end);
    #endif
    pst_data.shrinkJData(id);
  }


  inline void newJxn(jnid_t id, size_t max_size) {
    size_t tmp = jxn_data.newJData(max_size, true); 
    assert(tmp == id);
  }

  inline bool newUnion(jnid_t const id, size_t max_len, vid_t Xclude) {
    size_t sum = 0;
    std::vector<SortedRange> kid_itrs;
    kid_itrs.reserve(kids(id).size() + 1);
    for (auto itr = kids(id).cbegin(); itr != kids(id).cend(); itr++) {
      if (jxn(*itr).len != 0) {
        kid_itrs.emplace_back(jxn(*itr).begin(), jxn(*itr).end());
        sum += jxn(*itr).len;
      }
    }
    if (pst(id).len != 0) {
      kid_itrs.emplace_back(pst(id).begin(), pst(id).end());
      sum += pst(id).len;
    }

    max_len = std::min(max_len, sum);
    newJxn(id, max_len);
    bool success = heuristic_merge(jxn(id), max_len, kid_itrs, Xclude);
    if (success)
      jxn_data.shrinkJData(id);
    else
      jxn_data.deleteJData(id);
    return success;
  }

  inline JData<vid_t>& jxn(jnid_t id) { return jxn_data[id]; }
  inline JData<vid_t> const & jxn(jnid_t id) const { return jxn_data[id]; }

  inline size_t width(jnid_t id) const {
    return 1 + (id < jxn_data.size() ? jxn(id).size() : pst_weight(id));
  }
  

  /* JNODETABLE FAQS */
  inline void print(jnid_t id) const {
    printf("%6zu:w%6zu:pre%6zu:pst        ->[%4zu]\n",
      width(id), (size_t) pre_weight(id), (size_t) pst_weight(id),  (size_t) parent(id));
  }

  struct Facts {
    size_t vert_cnt;
    size_t edge_cnt;

    size_t width;
    size_t fill;

    size_t vert_height;
    size_t edge_height;
    size_t root_cnt;

    jnid_t halo_id;
    jnid_t core_id;

    Facts(JNodeTable const &jn);

    inline void print() const {
      printf("TREEFAQS: width:%zu\troots:%zu\n", width, root_cnt);
      printf("\tvheight:%zu\teheight:%zu\n", vert_height, edge_height);
      printf("\tverts:%zu\tedges:%zu\n", vert_cnt, edge_cnt);
      printf("\thalo:%zu\tcore:%zu\n", (size_t) halo_id, (size_t) core_id);
      printf("\tfill:%zu\n", fill);
    }
  };

  inline Facts getFacts() const { return Facts(*this); }
};

#include "jnode.cpp"

