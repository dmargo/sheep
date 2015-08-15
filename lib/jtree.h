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

#include <algorithm>
#include <fstream>
#include <vector>

#include "defs.h"
#include "graph_wrapper.h"
#include "jnode.h"

/* A JTREE represents the isomorphism between a graph and a chordal embedding (JNODES) via an INDEX.
 * In particular, JTree implements the algorithm to make a chordal embedding from a sequence isomorphism. */
class JTree {
private:
  std::vector<jnid_t> index; //Graph isomorphism; indexed by vid_t (vid_t to jnid_t mapping).
public:
  JNodeTable jnodes; //The chordal embedding; indexed (labeled) by jnid_t.

  inline jnid_t vid2jnid(vid_t X) const { return X < index.size() ? index[X] : INVALID_JNID; }

  inline size_t size() const { return jnodes.size(); }

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
      printf("%4zu:%-8zu", (size_t) id, (size_t) jnid2vid.at(id));
      jnodes.print(id);
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

#include "jtree.cpp"

