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

#include "jtree.h"

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

  if (opts.make_jxn && !jnodes.newUnion(current, opts.width_limit, X))
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

    /* XXX This code is hard to understand...it supports features that are likely to be cut.
     * If you are code-reading for the first time, you probably have no need to understand this. */
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

