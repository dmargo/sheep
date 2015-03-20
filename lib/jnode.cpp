#include "jnode.h"

JNodeTable::JNodeTable(jnid_t max_jnids, bool init_kids, size_t memory_limit) :
  nodes_state(State::ALLOCATED), end_id(0),
  max_id(max_jnids), nodes((JNode*)malloc(sizeof(JNode) * max_id)),
  kid_data(init_kids ? max_id : 0), pst_data(max_id, memory_limit), jxn_data(max_id, memory_limit),
  roots(max_id)
{
  if (nodes == nullptr)
    throw std::bad_alloc();
}

JNodeTable::JNodeTable(char const *filename, jnid_t max_jnids, bool init_kids, size_t memory_limit) :
  nodes_state(State::MAPPED), end_id(0), max_id(max_jnids), nodes(nullptr),
  kid_data(init_kids ? max_id : 0), pst_data(max_id, memory_limit), jxn_data(max_id, memory_limit),
  roots(max_id)
{
  int fd = open(filename, O_RDWR | O_CREAT | O_TRUNC, 00666);
  if (fd == -1)
    throw std::bad_alloc();

  if (posix_fallocate(fd, 0, sizeof(jnid_t) + sizeof(JNode) * max_id) != 0) {
    close(fd);
    throw std::bad_alloc();
  }

  char *const nodes_map = (char*)mmap(nullptr, sizeof(jnid_t) + sizeof(JNode) * max_id,
        PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (nodes_map == MAP_FAILED) {
    close(fd);
    throw std::bad_alloc();
  }
  close(fd);
  nodes = (JNode*)(nodes_map + sizeof(jnid_t));
}

JNodeTable::JNodeTable(char const *filename) :
  nodes_state(State::MAPPED), end_id(0), max_id(0), nodes(nullptr),
  kid_data(0), pst_data(0), jxn_data(0), roots(0)
{
  int fd = open(filename, O_RDWR);
  if (fd == -1)
    throw std::bad_alloc();

  struct stat buf;
  if (fstat(fd, &buf) == -1) {
    close(fd);
    throw std::bad_alloc();
  }
  max_id = (buf.st_size - sizeof(jnid_t)) / sizeof(JNode);

  char *const nodes_map = (char*)mmap(nullptr, sizeof(jnid_t) + sizeof(JNode) * max_id,
        PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (nodes_map == MAP_FAILED) {
    close(fd);
    throw std::bad_alloc();
  }
  close(fd);
  end_id = *((jnid_t*)nodes_map);
  nodes = (JNode*)(nodes_map + sizeof(jnid_t));

  kid_data = std::move(JDataTable<jnid_t>(max_id));
  makeKids();

  //XXX Currently there's nothing you can do with these, so why bother?
  #if 0
  pst_data = std::move(JDataTable<vid_t>(max_id, memory_limit));
  jxn_data = std::move(JDataTable<vid_t>(max_id, memory_limit));

  roots = UnionFind(max_id);
  for (jnid_t id = 0; id < size(); ++id)
    if (parent(id) != INVALID_JNID)
      roots.unify(id, parent(id));
  #endif
}

JNodeTable::JNodeTable(JNode *n, jnid_t end) :
  nodes_state(State::TEMPORARY), end_id(end), max_id(end_id), nodes(n),
  kid_data(max_id), pst_data(0), jxn_data(0), roots(0)
{
  if (nodes != nullptr)
    makeKids();
}

JNodeTable::JNodeTable(JNodeTable &&other) :
  nodes_state(other.nodes_state), end_id(other.end_id), max_id(other.max_id), nodes(other.nodes),
  kid_data(std::move(other.kid_data)), pst_data(std::move(other.pst_data)), jxn_data(std::move(other.jxn_data)),
  roots(std::move(other.roots))
{
  other.nodes_state = State::TEMPORARY;
  other.end_id = 0;
  other.max_id = 0;
  other.nodes = nullptr;
}

JNodeTable::JNodeTable(JNodeTable const &other) :
  nodes_state(State::ALLOCATED), end_id(other.end_id),
  max_id(other.max_id), nodes((JNode*)malloc(sizeof(JNode) * max_id)),
  kid_data(other.kid_data), pst_data(other.pst_data), jxn_data(other.jxn_data), roots(other.roots)
{
  if (nodes == nullptr)
    throw std::bad_alloc();
  std::memcpy(nodes, other.nodes, end_id);
}

JNodeTable::JNodeTable(JNodeTable const &other, jnid_t partial_end) :
  nodes_state(State::ALLOCATED), end_id(other.end_id),
  max_id(other.max_id), nodes((JNode*)malloc(sizeof(JNode) * max_id)),
  kid_data(other.kid_data, partial_end),
  pst_data(other.pst_data, partial_end),
  jxn_data(other.jxn_data, partial_end),
  roots(max_id)
{
  if (nodes == nullptr)
    throw std::bad_alloc();
  std::memcpy(nodes, other.nodes, partial_end);
  
  for (jnid_t id = 0; id < size(); ++id) {
    if (parent(id) < end_id)
      roots.unify(id, parent(id));
    else
      parent(id) = INVALID_JNID;
  }
}

JNodeTable::~JNodeTable() {
  if (nodes_state == State::ALLOCATED)
    free(nodes);
  else if (nodes_state == State::MAPPED) {
    char *const nodes_map = ((char*)nodes) - sizeof(jnid_t);
    *((jnid_t*)nodes_map) = end_id;
    munmap(nodes_map, sizeof(jnid_t) + sizeof(JNode*) * max_id);
  }
}

void JNodeTable::save(char const *const filename) {
  std::ofstream stream(filename, std::ios::binary | std::ios::trunc);
  stream.write((char*)&end_id, sizeof(jnid_t));
  stream.write((char*)nodes, sizeof(JNode) * max_id);
}

void JNodeTable::merge(JNodeTable const &lhs, JNodeTable const &rhs, bool const make_kids)
{
  assert(lhs.size() == rhs.size());

  for (jnid_t current = 0; current < lhs.size(); ++current) {
    jnid_t const tmp = newJNode();
    assert(tmp == current);

    if (make_kids) newKids(current, lhs.kids(current).size() + rhs.kids(current).size());

    //XXX For k-way merge, generalize to a JNodeTable list
    auto insert_kids = [current,make_kids,this](JNodeTable const &src)
    {
      for (jnid_t const kid : src.kids(current)) {
        if (!make_kids)
          adopt(kid, current);
        else
          meetKid(kid, current, src.pre_weight(kid));
      }
      pst_weight(current) += src.pst_weight(current);
    };
    insert_kids(lhs);
    insert_kids(rhs);

    if (make_kids)
      adoptKids(current);
  }
}

void mpi_merge_reduction(void *in, void *inout, int *len, MPI_Datatype *datatype) {
  JNodeTable lhs((JNodeTable::JNode*)in, *len);
  JNodeTable rhs((JNodeTable::JNode*)inout, *len);

  JNodeTable tmp(*len, false, 0);
  tmp.merge(lhs, rhs, false);

  assert(*len = tmp.end_id);
  memcpy(inout, tmp.nodes, sizeof(JNodeTable::JNode*) * tmp.end_id);
}

void JNodeTable::mpi_merge(bool const make_kids)
{
#ifndef USE_PRE_WEIGHT
  int count = 2;
  int blocklengths[2] = { 1, 1 };
  MPI_Aint offsets[2] = { offsetof(JNode, parent), offsetof(JNode, pst_weight) };
  MPI_Datatype types[2] = { MPI_UNSIGNED, MPI_UNSIGNED };
#else
  int count = 3;
  int blocklengths[3] = { 1, 1, 1 };
  MPI_Aint offsets[3] = { offsetof(JNode, parent), offsetof(JNode, pst_weight), offsetof(JNode, pre_weight) };
  MPI_Datatype types[3] = { MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED };
#endif
  MPI_Datatype mpi_jnode_type;
  MPI_Type_create_struct(count, blocklengths, offsets, types, &mpi_jnode_type);
  MPI_Type_commit(&mpi_jnode_type);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  JNode *outbuf = rank != 0 ? nullptr : (JNode*)malloc(sizeof(JNode) * max_id);

  MPI_Op reduce_op;
  MPI_Op_create(mpi_merge_reduction, 0, &reduce_op);
  MPI_Reduce(nodes, outbuf, end_id, mpi_jnode_type, reduce_op, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (nodes_state == State::ALLOCATED)
      std::swap(nodes, outbuf);
    else
      memcpy(nodes, outbuf, sizeof(JNode) * end_id);
    free(outbuf);
  }
}



/*XXX This stuff doesn't really belong here and needs to be moved. */
struct SortedRange {
  vid_t const *itr;
  vid_t const *end;

  inline SortedRange(vid_t const *begin, vid_t const *e) : itr(begin), end(e) {}
  inline bool operator>(SortedRange const &rhs) const { return (*itr > *rhs.itr); }
  inline size_t size() const { return std::distance(itr,end); }
};

bool balance_line_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
{
  for (auto end = kid_itrs.end(); end != kid_itrs.begin();) {
    // Find range of minimums (linear search)
    vid_t min = std::numeric_limits<vid_t>::max();
    std::vector<SortedRange>::iterator frst, last;
    for (auto cur = kid_itrs.begin(); cur != end; cur++) {
      if (*cur->itr < min) {
        min = *cur->itr;
        frst = last = cur;
      } else if (*cur->itr == min) {
        last = cur;
      }
    }

    // Insert min; always need to exclude one elimination element.
    if (min != Xclude) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(min);
    }
    
    // Increment min (linear search)
    for (; frst != last + 1 && frst != end; frst++)
      if (*frst->itr == min)
        if (++(frst->itr) == frst->end)
          std::iter_swap(--end, frst--);
  }
  return true;
}

/* heap merge may outscale balance-line merge for something like kid_itrs.size() > 32.
 * The likelihood of this case increases with graph density.
 * XXX: TBD whether this is significant for Big graphs.
 * XXX: I've seen at least one case where it made a HUGE difference. Needs to be benchmarked.
 * does priority_queue virtualize std::greater calls? */
bool heap_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
{
  std::priority_queue< SortedRange, std::vector<SortedRange>, std::greater<SortedRange> >
      q(std::greater<SortedRange>(), std::move(kid_itrs));
  vid_t prev = INVALID_VID;

  while (!q.empty()) {
    SortedRange cur = q.top();
    q.pop();

    vid_t const min = *cur.itr++;
    if (min != Xclude && min != prev) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(min);
      prev = min;
    }

    if (cur.itr != cur.end)
      q.push(cur);
  }
  return true;
}

bool asymmetric_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
{
  assert(kid_itrs.size() == 2);
  SortedRange big = kid_itrs[0];
  SortedRange small = kid_itrs[1];
  if (big.size() < small.size())
    std::swap(big,small);
  if (big.size() < small.size() * 8)
    return balance_line_merge(new_data, kid_itrs, Xclude, max_len);

  for (; small.itr != small.end; small.itr++) {
    auto *big_middle = std::lower_bound(big.itr, big.end, *small.itr);

    for (; big.itr != big_middle; big.itr++) {
      if (*big.itr != Xclude) {
        if (new_data.len + 1 > max_len)
          return false;
        new_data.push_back(*big.itr);
      }
    }

    if (*small.itr != Xclude && *small.itr != *big.itr) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(*small.itr);
    }
  }

  for (; big.itr != big.end; big.itr++) {
    if (*big.itr != Xclude) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(*big.itr);
    }
  }

  return true;
}

bool JNodeTable::newUnion(jnid_t const id, vid_t Xclude, size_t max_len) {
  std::vector<SortedRange> kid_itrs;
  kid_itrs.reserve(kids(id).size() + 1);
  size_t sum = 0;

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

  bool success;
  if (kid_itrs.size() < 32)
    success = balance_line_merge(jxn(id), kid_itrs, Xclude, max_len);
  else
    success = heap_merge(jxn(id), kid_itrs, Xclude, max_len);

  if (success)
    jxn_data.shrinkJData(id);
  else
    jxn_data.deleteJData(id);

  return success;
}
