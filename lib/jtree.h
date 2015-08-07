#pragma once

#include <algorithm>
#include <fstream>
#include <vector>

#include "graph_wrapper.h"
#include "jnode.h"
#include "stdafx.h"

/* A JTREE represents the isomorphism between a graph and a chordal embedding (JNODES) via an INDEX.
 * In particular, JTree implements the algorithm to make a chordal embedding from a sequence isomorphism. */
class JTree {
private:
  std::vector<jnid_t> index; //Graph isomorphism; indexed by vid_t (vid_t to jnid_t mapping).
public:
  JNodeTable jnodes; //The chordal embedding; indexed (labeled) by jnid_t.

  inline size_t size() const { return jnodes.size(); }
  
  inline jnid_t vid2jnid(vid_t X) const { return X < index.size() ? index[X] : INVALID_JNID; }

  inline std::vector<vid_t> get_sequence() const {
    std::vector<vid_t> seq(size());
    for (vid_t X = 0; X != index.size(); ++X)
      if (index[X] != INVALID_JNID)
        seq.at(index[X]) = X;
    return seq;
  }

  inline void print() const {
    std::vector<vid_t> jnid2vid = get_sequence();
    for (jnid_t id = 0; id != size(); ++id) {
      printf("%4u:%-8u", id, jnid2vid.at(id));
      jnodes.print(id);
    }
  }

  template <typename GraphType>
  inline void write_isomorphism(GraphType const &graph, char const *const filename) const {
    std::ofstream stream(filename);

    std::vector<vid_t> jnid2vid = get_sequence();
    for (jnid_t X_id = 0; X_id != size(); ++X_id) {
      for (auto eitr = graph.getEdgeItr(jnid2vid.at(X_id)); !eitr.isEnd(); ++eitr) {
        jnid_t const Y_id = vid2jnid(*eitr);
        assert(Y_id != INVALID_JNID);
        stream << X_id << ' ' << Y_id << std::endl;
      }
    }
  }


  /* JTREE OPTIONS */
  struct Options {
    bool verbose;

    bool make_pad;  // make jnodes for verts with no edges
    bool make_kids; // make out-tree (child) pointers
    bool make_pst;  // make post-neighbor edge table
    bool make_jxn;  // make fill-neighbor edge table

    size_t memory_limit;  // limit the maximum memory used for pst and jxn tables
    size_t width_limit;   // defer vertices of width > width_limit to the end of the input sequence
    bool find_max_width;  // quit when we find the max width (treewidth) of the sequence

    bool do_rooting;
    size_t rooting_limit;

    Options() :
      verbose(false),
      make_pad(true), make_kids(false), make_pst(false), make_jxn(false),
      memory_limit(1 * GIGA), width_limit((size_t)-1), find_max_width(false),
      do_rooting(false), rooting_limit(0) {}

    bool isDefault() const {
      return
        verbose == false &&
        make_pad == true && make_kids == false && make_pst == false && make_jxn == false &&
        memory_limit == 1 * GIGA && width_limit == 0 && find_max_width == false &&
        do_rooting == false && rooting_limit == 0;
    }

    bool isValid() const {
      return
        (make_jxn ? (make_kids && make_pst) : true) &&
        (width_limit != (size_t)-1 ? make_jxn : true) &&
        (find_max_width ? make_jxn : true) &&
        (do_rooting ? make_jxn : true) &&
        (rooting_limit != 0 ? do_rooting : true);
    }
  };

  /* CONSTRUCTORS */
  template <typename GraphType>
  inline JTree(GraphType const &graph, std::vector<vid_t> const &seq,
    Options opts = Options()) : index(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID),
    jnodes(opts.make_pad ? seq.size() : graph.getNodes(),
           opts.make_kids,
           opts.make_pst || opts.make_jxn ? opts.memory_limit : 0)
  {
    if (opts.isDefault())
      insertSequence(graph, seq);
    else
      insertSequence(graph, seq, opts);
  }

  template <typename GraphType>
  inline JTree(GraphType const &graph, std::vector<vid_t> const &seq, char const *const filename,
    Options opts = Options()) : index(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID),
    jnodes(filename,
           opts.make_pad ? seq.size() : graph.getNodes(),
           opts.make_kids,
           opts.make_pst ? opts.memory_limit : 0)
  {
    if (opts.isDefault())
      insertSequence(graph, seq);
    else
      insertSequence(graph, seq, opts);
  }

  /* open constructor */
  inline JTree(std::vector<vid_t> const &seq, char const *const filename) :
    index(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID), jnodes(filename)
  {
    for (jnid_t id = 0; id != seq.size(); ++id)
      index[seq[id]] = id;
  }

  /* copy constructors */
  inline JTree(JTree &&other) = default;

  inline JTree(JTree const &other, jnid_t partial_end) :
    index(other.index), jnodes(other.jnodes, partial_end)
  {
    std::replace_if(index.begin(), index.end(),
      [partial_end](jnid_t const id) { return id < partial_end; }, INVALID_JNID);
  }

  /* special members */
  ~JTree() = default;

  JTree() = delete;
  JTree(JTree const &other) = delete;
  JTree& operator=(JTree &&other) = delete;
  JTree& operator=(JTree const &other) = delete;

private:
  inline void insert(vid_t X, jnid_t id) {
    assert(index.at(X) == INVALID_JNID);
    index.at(X) = id;
  }

  template <typename GraphType>
  jnid_t insert(GraphType const &graph, vid_t X);

  template <typename GraphType>
  void insertSequence(GraphType const &graph, std::vector<vid_t> const &seq);

  template <typename GraphType>
  jnid_t insert(GraphType const &graph, vid_t X, Options opts);

  template <typename GraphType>
  void insertSequence(GraphType const &graph, std::vector<vid_t> const &seq, Options opts);

public:
  template <typename GraphType>
  bool isValid(GraphType const &graph, std::vector<vid_t> const &seq, Options opts = Options()) const;
};





//XXX DRY, but these non-parameterized versions make like a 10% performance difference.
//They are also much easier to read and understand, so they serve a documentary purpose.
template <typename GraphType>
jnid_t JTree::insert(GraphType const &graph, vid_t const X)
{
  jnid_t const current = jnodes.newJNode();

  if (graph.isNode(X)) {
    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const nbr = *eitr;
      jnid_t const nbr_id = index.at(nbr);

      // PREORDER edge
      if (nbr_id != INVALID_JNID)
        jnodes.adopt(nbr_id, current);
      // POSTORDER edge
      else if (nbr != X)
        ++jnodes.pst_weight(current);
    }
  }

  insert(X, current);
  return current;
}

template <typename GraphType>
void JTree::insertSequence(GraphType const &graph, std::vector<vid_t> const &seq)
{
  for (vid_t const X : seq)
    insert(graph, X);
}

// Parameterized insert
template <typename GraphType>
jnid_t JTree::insert(GraphType const &graph, vid_t const X, Options const opts)
{
  jnid_t const current = jnodes.newJNode();
  if (opts.make_kids) jnodes.newKids(current, graph.isNode(X) ? graph.getDeg(X) : 0);
  if (opts.make_pst) jnodes.newPst(current, graph.isNode(X) ? graph.getDeg(X) : 0);

  if (graph.isNode(X)) {
    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const nbr = *eitr;
      vid_t const nbr_id = index.at(nbr);

      // PREORDER edge
      if (nbr_id != INVALID_JNID) {
        if (!opts.make_kids)
          jnodes.adopt(nbr_id, current);
        else
          jnodes.meetKid(nbr_id, current, 1);  
      }
      // POSTORDER edge
      else if (nbr != X) {
        if (++jnodes.pst_weight(current) > opts.width_limit)
          goto FAILURE;
        else if (opts.make_pst)
          jnodes.pst(current).push_back(nbr);
      }
    }
  }

  if (opts.make_pst)
    jnodes.cleanPst(current);

  if (opts.make_jxn && !jnodes.newUnion(current, X, opts.width_limit))
    goto FAILURE;

  //XXX This cannot be revoked, so it must be deferred until now.
  if (opts.make_kids)
    jnodes.adoptKids(current);
  
  insert(X, current);
  return current;
  
FAILURE:
  jnodes.deleteJNode(current);
  return INVALID_JNID;
}

template <typename GraphType>
void JTree::insertSequence(GraphType const &graph, std::vector<vid_t> const &seq, Options const opts) {
  assert(opts.isValid());
  if (opts.verbose) printf("Constructing JTree.");

  auto seq_itr = seq.cbegin();
  std::vector<vid_t> wide_seq;
  for (size_t current_width = 0; seq_itr != seq.cend(); ++seq_itr)
  {
    if (opts.verbose && size() % 1000 == 0) {
      size() % 1000000 == 0 ? printf("%zu", size() / 1000000) : printf(".");
      fflush(stdout);
    }

    vid_t const X = *seq_itr;
    if (!opts.make_pad && !graph.isNode(X)) continue;
    jnid_t const current = insert(graph, X, opts);

    if (opts.find_max_width) {
      if (current != INVALID_JNID)
        current_width = std::max(current_width, jnodes.width(current));
      if (current_width >= wide_seq.size() + std::distance(seq_itr, seq.cend()))
        return;
    }

    if (current == INVALID_JNID)
      wide_seq.push_back(X);
    else if (opts.do_rooting && jnodes.width(current) == wide_seq.size() + std::distance(seq_itr, seq.cend())) {
      ++seq_itr;
      break;
    }
  }

  // DRY; but this code is very performant in the special case it covers.
  // What I suspect is that the case will gradually become more special until it is distinguishable.
  // At the very least it should be extracted into helper functions, which will be much easier
  // if I can find a good way to reinsert wide_seq into seq.
  // XXX pre_weight is currently broken for this case.
  auto wide_itr = wide_seq.cbegin();
  if (wide_itr != wide_seq.cend() || seq_itr != seq.cend()) {
    size_t remaining = wide_seq.size() + std::distance(seq_itr, seq.cend()) - 1;

    vid_t X = wide_itr != wide_seq.cend() ? *wide_itr++ : *seq_itr++;
    jnid_t current = jnodes.newJNode();

    jnodes.newKids(current, jnodes.size());
    for (jnid_t kid = 0; kid < jnodes.size(); ++kid) {
      if (jnodes.parent(kid) == INVALID_JNID && kid != current) {
        jnodes.parent(kid) = current;
        jnodes.kids(current).push_back(kid);
      }
    }
    jnodes.adoptKids(current);

    jnodes.newPst(current, graph.isNode(X) ? graph.getDeg(X) : 0);
    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const nbr = *eitr;
      if (index.at(nbr) == INVALID_JNID && nbr != X) {
        ++jnodes.pst_weight(current);
        jnodes.pst(current).push_back(nbr);
      }
    }
    jnodes.adoptKids(current);

    // The special case is that jxn == remaining, so this (very large) union is trivial
    jnodes.newJxn(current, remaining);
    for (auto inner_wide_itr = wide_itr; inner_wide_itr != wide_seq.cend(); ++inner_wide_itr)
      jnodes.jxn(current).push_back(*inner_wide_itr);
    for (auto inner_seq_itr = seq_itr; inner_seq_itr != seq.cend(); ++inner_seq_itr)
      jnodes.jxn(current).push_back(*inner_seq_itr);
    std::sort(jnodes.jxn(current).begin(), jnodes.jxn(current).end());

    // Done; insert
    insert(X, current);
    if (opts.find_max_width) return;

    // Once rooted, subsequent vertices are also trivial
    while (wide_itr != wide_seq.cend() || seq_itr != seq.cend()) {
      remaining--;

      vid_t X = wide_itr != wide_seq.cend() ? *wide_itr++ : *seq_itr++;
      jnid_t const previous = current;
      current = jnodes.newJNode();

      jnodes.newKids(current, 1);
      jnodes.parent(previous) = current;
      jnodes.kids(current).push_back(previous);
      jnodes.adoptKids(current);

      jnodes.newPst(current, graph.isNode(X) ? graph.getDeg(X) : 0);
      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const nbr = *eitr;
        if (index.at(nbr) == INVALID_JNID && nbr != X) {
          ++jnodes.pst_weight(current);
          jnodes.pst(current).push_back(nbr);
        }
      }
      jnodes.cleanPst(current);

      jnodes.newJxn(current, remaining);
      for (auto itr = jnodes.jxn(previous).begin(); itr != jnodes.jxn(previous).end(); ++itr)
        if (*itr != X)
          jnodes.jxn(current).push_back(*itr);

      insert(X, current);
    }

    assert(remaining == 0);
  }

  if (opts.verbose) {
    printf("done\n");
    if (size() < graph.getNodes())
      printf("WARNING insertSequence(): seq is not a total order on the graph's remaining vertices.\n");
    else if (size() > graph.getNodes())
      printf("WARNING insertSequence(): seq adds 0-degree vertices not in reference graph.\n");
  }
}

#ifndef NDEBUG
#define FAIL_IF(bool_exp) assert(!(bool_exp))
#else
#define FAIL_IF(bool_exp) if (bool_exp) return false
#endif
template <typename GraphType>
bool JTree::isValid(GraphType const &graph, std::vector<vid_t> const &seq, Options const opts) const{
  JNodeTable::Facts const faq = jnodes.getFacts();

  jnid_t valid_indices = 0;
  for (jnid_t id : index)
    if (id != INVALID_JNID)
      ++valid_indices;
  FAIL_IF(valid_indices != faq.vert_cnt);
  
#ifndef NDEBUG
#pragma omp parallel for
#endif
  for (size_t i = 0; i < seq.size(); ++i) {
    vid_t const X = seq.at(i);
    if (!opts.make_pad && !graph.isNode(X)) continue;

    jnid_t current = index.at(X);
    FAIL_IF(current == INVALID_JNID);
    FAIL_IF(current < 0);
    FAIL_IF(current >= jnodes.size());

    FAIL_IF(opts.make_pst && jnodes.pst(current).binary_search(X)); //no self-edges
    FAIL_IF(opts.make_jxn && jnodes.jxn(current).binary_search(X)); //common union bug
    if (opts.make_kids)
      for (jnid_t kid : jnodes.kids(current))
        FAIL_IF(jnodes.parent(kid) != current);
    
    if (graph.isNode(X)) {
      for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
        vid_t const nbr = *eitr;
        jnid_t nbr_id = index.at(X);

        if (nbr_id < current) {
          for (size_t step = 0; nbr_id != current; step++) {
            FAIL_IF(step > faq.vert_height);

            FAIL_IF(nbr_id == INVALID_JNID);
            FAIL_IF(nbr_id < 0);
            FAIL_IF(nbr_id >= jnodes.size());
            FAIL_IF(opts.make_jxn && !jnodes.jxn(nbr_id).binary_search(X)); //IMPORTANT

            nbr_id = jnodes.parent(nbr_id);
          }
        }
        else if (nbr_id > current) {
          FAIL_IF(opts.make_pst && !jnodes.pst(current).binary_search(nbr));
          FAIL_IF(opts.make_jxn && !jnodes.jxn(current).binary_search(nbr));
        }
      }
    }

    for (size_t step = 0; jnodes.parent(current) != INVALID_JNID; step++) {
      FAIL_IF(step > faq.vert_height);
      
      FAIL_IF(current < 0);
      FAIL_IF(current >= jnodes.size());

      current = jnodes.parent(current);
    }
  }
  return true;
}
#undef FAIL_IF

