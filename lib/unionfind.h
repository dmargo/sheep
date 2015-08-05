#include <numeric>
#include <vector>

// T should be a type that can index an array.
template <typename T>
class FastUnionFind {
protected:
  std::vector<T> parent;
  std::vector<char> rank;
  //XXX The max value I've seen for rank is 8 on very large graphs.
  //I believe it needs to kept separate from parent for alignment savings, though this is not? cache-optimal.
  //It is also very sparse. However, it saves measurable time. Not sure what to trade-off.

  //XXX This is currently the innermost loop of the most basic algorithm.
  //Like 9% of runtime is spent in std::vector::operator[]
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

