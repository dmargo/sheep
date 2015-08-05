#pragma once
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>
#include <vector>

#include "jdata.h"
#include "unionfind.h"

#include <fcntl.h>
#include <mpi.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

typedef vid_t jnid_t;
#define INVALID_JNID ((jnid_t)-1)

typedef FastUnionFind<jnid_t> UnionFind;

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

    for (jnid_t const kid : kids(id)) {
      parent(kid) = INVALID_JNID;
      pre_weight(kid) = 0;
    }
    kid_data.deleteJData(id);
    pst_data.deleteJData(id);
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

  //XXX Complicated; consider exporting this to its own merge.h
  bool newUnion(jnid_t id, vid_t Xclude, size_t max_len = (size_t)-1);

  inline JData<vid_t>& jxn(jnid_t id) { return jxn_data[id]; }
  inline JData<vid_t> const & jxn(jnid_t id) const { return jxn_data[id]; }

  inline size_t width(jnid_t id) const {
    return 1 + (id < jxn_data.size() ? jxn(id).size() : pst_weight(id));
  }
  


  /* JNODETABLE FAQS */
  inline void print(jnid_t id) const {
    printf("%6zu:w%6u:pre%6u:pst        ->[%4u]\n",
      width(id), pre_weight(id), pst_weight(id),  parent(id));
  }

  struct Facts {
    size_t vert_cnt;
    size_t edge_cnt;

    size_t width;
    long long unsigned fill;

    long long unsigned vert_height;
    long long unsigned edge_height;
    size_t root_cnt;

    jnid_t halo_id;
    jnid_t core_id;

    std::vector<size_t> distr;

    Facts(JNodeTable const &jn);

    inline void print() const {
      printf("TREEFAQS: width:%zu\troots:%zu\n", width, root_cnt);
      printf("\tvheight:%llu\teheight:%llu\n", vert_height, edge_height);
      printf("\tverts:%zu\tedges:%zu\n", vert_cnt, edge_cnt);
      printf("\thalo:%u\tcore:%u\n", halo_id, core_id);
      printf("\tfill:%llu\n", fill);
    }

    inline void print_distribution() const {
      for (size_t i = 0; i != distr.size(); ++i)
        if (distr.at(i) != 0)
          printf("%zu : %zu\n", i, distr.at(i));
    }
  };

  inline Facts getFacts() const { return Facts(*this); }
};

inline JNodeTable::Facts::Facts(JNodeTable const &jnodes) :
  vert_cnt(0), edge_cnt(0), width(0), fill(0),
  vert_height(0), edge_height(0), root_cnt(0),
  halo_id(INVALID_JNID), core_id(INVALID_JNID), 
  distr()
{
  std::vector<long long unsigned> vheight(jnodes.size(), 0);
  std::vector<long long unsigned> eheight(jnodes.size(), 0);

  // Ascending pass; it is natural to compute most facts here.
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    jnid_t const par_id = jnodes.parent(id);

    vert_cnt++;
    edge_cnt += jnodes.pst_weight(id);
    width = std::max(width, jnodes.width(id));
    fill += jnodes.width(id) - jnodes.pst_weight(id) - 1;

    vheight.at(id)++;
    eheight.at(id) += jnodes.pst_weight(id);
    if (par_id != INVALID_JNID) {
      vheight.at(par_id) = std::max<long long unsigned>(vheight.at(par_id), vheight.at(id));
      eheight.at(par_id) = std::max<long long unsigned>(eheight.at(par_id), eheight.at(id));
    }
    else {
      vert_height = std::max<long long unsigned>(vert_height, vheight.at(id));
      edge_height = std::max<long long unsigned>(edge_height, eheight.at(id));
      root_cnt++;
    }

    if (halo_id == INVALID_JNID && jnodes.width(id) > 3)
      halo_id = id;
    if (core_id == INVALID_JNID && jnodes.width(id) >= width)
      core_id = id;

#if 0
    if (distr.size() < jnodes.kids(id).size() + 1)
      distr.resize(jnodes.kids(id).size() + 1);
    distr.at(jnodes.kids(id).size()) += 1;
#endif
  }

#if 0
  // Guaranteed to clear and free vectors.
  std::vector<long long unsigned>().swap(vheight);
  std::vector<long long unsigned>().swap(eheight);

  std::vector<long long unsigned> vdepth(jnodes.size(), 0);
  long long unsigned trivial_cost = edge_cnt * ceil(log2(vert_cnt));
  long long unsigned linear_cost = 0;
  long long unsigned tree_cost = vert_cnt * ceil(log2(vert_cnt));

  // Descending pass; generally for depth-oriented facts.
  for(jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
    jnid_t const par_id = jnodes.parent(id);
    if (par_id != INVALID_JNID)
      vdepth.at(id) = vdepth.at(par_id) + 1;
 
    if (jnodes.size() - 1 - id != 0)
      linear_cost += jnodes.pst_weight(id) * ceil(log2(jnodes.size() - 1 - id));
    if (vdepth.at(id) != 0)
      tree_cost += jnodes.pst_weight(id) * ceil(log2(vdepth.at(id)));

    if (distr.size() < vdepth.at(id) + 1)
      distr.resize(vdepth.at(id) + 1);
    distr.at(vdepth.at(id)) += 1;
  }

  //TODO: Move this into official print, assuming you're not discontinuing.
  printf("\nCOSTS:\n");
  printf("\ttrivial: %llu\n", trivial_cost);
  printf("\tlinear: %llu\n", linear_cost);
  printf("\ttree: %llu\n", tree_cost);
#endif
}
