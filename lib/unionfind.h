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

#include <numeric>
#include <vector>

// T should be a type that can index an array.
template <typename T>
class FastUnionFind {
protected:
  std::vector<T> parent;
  std::vector<char> rank;
  //XXX If rank overflows it is NOT a correctness issue. It only affects performance.
  //The maximum value I've seen for rank is 8 on very large graphs.
  //I believe it needs to kept separate from parent for aligned space saving.
  //It is also very sparse. However, it improves performance significantly.

  //XXX This is currently the innermost loop of the most basic algorithm.
  //Like 9% of runtime is spent here in std::vector::operator[]
  inline T find_root(T const element) {
    assert (element < parent.size());

    T itr = element;
    while (rank[itr] < rank[parent[itr]])
      itr = parent[itr];
    T const root = itr;

    // Path compression
    itr = element;
    while (itr != root) {
      T const next = parent[itr];
      parent[itr] = root;
      itr = next;
    }

    return root;
  }

public:
  inline FastUnionFind(T const universe) : parent(universe), rank(universe) {
    std::iota(parent.begin(), parent.end(), (T)0);
  }

  inline FastUnionFind(FastUnionFind const &other, T partial_end) :
    parent(other.parent.size()), rank(other.rank.size())
  {
    assert(partial_end <= other.parent.size());
    std::copy(other.parent.cbegin(), other.parent.cend() + partial_end, parent.begin());
    std::iota(parent.begin() + partial_end, parent.end(), partial_end);
  }

  inline T find(T const element) {
    return parent[find_root(element)];
  }

  inline T unify(T const lesser, T const greater) {
    assert(lesser < greater);

    T const greater_root = find_root(greater);
    T const lesser_root = find_root(lesser);
    T const old_parent = parent[lesser_root];

    if (lesser_root != greater_root) {
      if (rank[lesser_root] > rank[greater_root]) {
        parent[lesser_root] = greater;
        parent[greater_root] = lesser_root;
      }
      else {
        assert(parent[greater_root] == greater);
        parent[lesser_root] = greater_root;
        if (rank[lesser_root] == rank[greater_root])
          rank[greater_root] += 1;
      }
    }
    return old_parent;
  }
};



//XXX: SimpleUnionFind uses marginally less memory than FastUnionFind.
template <typename T>
class SimpleUnionFind {
protected:
  std::vector<T> membership;

public:
  inline SimpleUnionFind() : membership() {}

  inline SimpleUnionFind(T universe) : membership(universe) {
    std::iota(membership.begin(), membership.end(), (T)0);
  }

  inline SimpleUnionFind(SimpleUnionFind const &other, T partial_end) :
    membership(other.membership.size())
  {
    assert(partial_end <= membership.size());
    std::copy(other.membership.cbegin(), other.membership.cbegin() + partial_end, membership.begin());
    std::iota(membership.begin() + partial_end, membership.end(), partial_end);
  }

  inline T find(T const element) {
    assert(element < membership.size());

    T itr = element;
    while (itr != membership[itr])
        itr = membership[itr];
    T const root = itr;

    // Path compression
    itr = element;
    while (itr != root) {
      T const next = membership[itr];
      membership[itr] = root;
      itr = next;
    }

    return root;
  }

  inline T unify(T const child, T const parent) {
    assert(child < membership.size());
    assert(parent < membership.size());

    T old_parent = membership[child];
    membership[child] = parent;
    return old_parent;
  }

  inline void revoke(T const child) {
    assert(child < membership.size());
    membership[child] = child;
  }
};

