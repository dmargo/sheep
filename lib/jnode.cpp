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

#include "jnode.h"

#include <cstring>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

/* 
 * CONSTRUCTORS AND I/O
 */
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

  makeKids();
}

JNodeTable::JNodeTable(JNode *n, jnid_t end) :
  nodes_state(State::TEMPORARY), end_id(end), max_id(end_id), nodes(n),
  kid_data(0), pst_data(0), jxn_data(0), roots(0)
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


/*
 * TREE MERGING METHODS
 */
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

template <bool make_kids>
void mpi_merge_reduction(void *in, void *inout, int *len, MPI_Datatype *datatype) {
  JNodeTable lhs((JNodeTable::JNode*)in, *len);
  JNodeTable rhs((JNodeTable::JNode*)inout, *len);

  JNodeTable tmp(*len, make_kids, 0);
  tmp.merge(lhs, rhs, make_kids);
  memcpy(inout, tmp.nodes, sizeof(JNodeTable::JNode) * tmp.end_id);
}

void JNodeTable::mpi_merge(bool const make_kids)
{
  MPI_Datatype MPI_jnid_t = sizeof(jnid_t) == 4 ? MPI_UINT32_T : MPI_UINT64_T;
  MPI_Datatype MPI_esize_t = sizeof(esize_t) == 4 ? MPI_UINT32_T : MPI_UINT64_T;
#ifndef USE_PRE_WEIGHT
  int count = 2;
  int blocklengths[2] = { 1, 1 };
  MPI_Aint offsets[2] = { offsetof(JNode, parent), offsetof(JNode, pst_weight) };
  MPI_Datatype types[2] = { MPI_jnid_t, MPI_esize_t };
#else
  int count = 3;
  int blocklengths[3] = { 1, 1, 1 };
  MPI_Aint offsets[3] = { offsetof(JNode, parent), offsetof(JNode, pst_weight), offsetof(JNode, pre_weight) };
  MPI_Datatype types[3] = { MPI_jnid_t, MPI_esize_t, MPI_esize_t };
#endif
  MPI_Datatype mpi_jnode_type;
  MPI_Type_create_struct(count, blocklengths, offsets, types, &mpi_jnode_type);
  MPI_Type_commit(&mpi_jnode_type);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  JNode *outbuf = rank != 0 ? nullptr : (JNode*)malloc(sizeof(JNode) * max_id);

  MPI_Op reduce_op;
  if (make_kids)
    MPI_Op_create(mpi_merge_reduction<true>, 0, &reduce_op);
  else
    MPI_Op_create(mpi_merge_reduction<false>, 0, &reduce_op);
  MPI_Reduce(nodes, outbuf, end_id, mpi_jnode_type, reduce_op, 0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    if (nodes_state == State::ALLOCATED)
      std::swap(nodes, outbuf);
    else
      memcpy(nodes, outbuf, sizeof(JNode) * end_id);
    free(outbuf);
  }
}


/* 
 * FAQ
 */
inline JNodeTable::Facts::Facts(JNodeTable const &jnodes) :
  vert_cnt(0), edge_cnt(0), width(0), fill(0),
  vert_height(0), edge_height(0), root_cnt(0),
  halo_id(INVALID_JNID), core_id(INVALID_JNID) 
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
  }
}

