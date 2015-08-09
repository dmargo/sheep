#pragma once

#include <algorithm>
#include <fstream>
#include <vector>

#include <mpi.h>
#include <parallel/algorithm>

#include "defs.h"
#include "readerwriter.h"

/* SEQUENCE CONSTRUCTORS */
template <typename GraphType>
std::vector<vid_t> defaultSequence(GraphType const &graph) {
  std::vector<vid_t> seq;
  seq.reserve(graph.getNodes());
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr)
    seq.push_back(*nitr);
  return seq;
}

template <typename GraphType>
std::vector<vid_t> degreeSequence(GraphType const &graph) {
  std::vector<vid_t> seq = defaultSequence(graph);
  __gnu_parallel::sort(seq.begin(), seq.end(), [&graph](vid_t const lhs, vid_t const rhs)
  {
    if (graph.getDeg(lhs) != graph.getDeg(rhs))
      return graph.getDeg(lhs) < graph.getDeg(rhs);
    else
      return lhs < rhs;
  });
  return seq;
}

template <typename GraphType>
std::vector<vid_t> mpiSequence(GraphType const &graph) {
  MPI_Datatype MPI_vid_t = sizeof(vid_t) == 4 ? MPI_UINT32_T : MPI_UINT64_T;
  MPI_Datatype MPI_esize_t = sizeof(esize_t) == 4 ? MPI_UINT32_T : MPI_UINT64_T;

  vid_t max_vid = 0;
  vid_t local_max = graph.getMaxVid();
  MPI_Allreduce((void*)&local_max, (void*)&max_vid, 1, MPI_vid_t, MPI_MAX, MPI_COMM_WORLD);

  esize_t *degree = (esize_t*)malloc((max_vid + 1) * sizeof(esize_t));
  esize_t *local_degree = (esize_t*)calloc(max_vid + 1, sizeof(esize_t));
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr)
    local_degree[*nitr] = graph.getDeg(*nitr);
  MPI_Allreduce((void*)local_degree, (void*)degree, max_vid, MPI_esize_t, MPI_SUM, MPI_COMM_WORLD);
  free(local_degree);
  
  std::vector<vid_t> seq;
  for (vid_t X = 0; X != max_vid + 1; ++X)
    if (degree[X] != 0)
      seq.push_back(X);

  __gnu_parallel::sort(seq.begin(), seq.end(), [degree](vid_t const lhs, vid_t const rhs)
  {
    if (degree[lhs] != degree[rhs])
      return degree[lhs] < degree[rhs];
    else
      return lhs < rhs;
  });
  free(degree);
  return seq;
}

std::vector<vid_t> fileSequence(char const *const filename) {
  //XXX DRY; fix this later
  std::vector<vid_t> degree;
  if (strcmp(".dat", filename + strlen(filename) - 4) == 0) {
    XS1Reader reader(filename);

    vid_t X,Y;
    while(reader.read(X,Y)) {
      size_t const required_size = std::max(X,Y) + 1;
      if (degree.size() < required_size) degree.resize(required_size);
      degree[X] += 1;
      degree[Y] += 1;
    }
  } else if (strcmp(".net", filename + strlen(filename) - 4) == 0) {
    SNAPReader reader(filename);

    vid_t X,Y;
    while(reader.read(X,Y)) {
      size_t const required_size = std::max(X,Y) + 1;
      if (degree.size() < required_size) degree.resize(required_size);
      degree[X] += 1;
      degree[Y] += 1;
    }

  } else {
    printf("sequence.h:netSequence(): Unsupported file type.\n");
    exit(1);
  }

  std::vector<vid_t> seq;
  for (vid_t X = 0; X != degree.size(); ++X)
    if (degree[X] != 0)
      seq.push_back(X);

  __gnu_parallel::sort(seq.begin(), seq.end(), [&degree](vid_t const lhs, vid_t const rhs)
  {
    if (degree[lhs] != degree[rhs])
      return degree[lhs] < degree[rhs];
    else
      return lhs < rhs;
  });
  return seq;
}



/* SEQUENCE I/O */
void writeBinarySequence(std::vector<vid_t> const &seq, char const *const filename) {
  std::ofstream stream(filename, std::ios::binary | std::ios::trunc);
  
  size_t const size = seq.size();
  stream.write((char*)(&size), sizeof(size_t));

  stream.write((char*)seq.data(), seq.size() * sizeof(vid_t));
}

std::vector<vid_t>  readBinarySequence(char const *const filename) {
  std::ifstream stream(filename, std::ios::binary);
  //XXX Test if file exists.

  size_t size = 0;
  stream.read((char*)&size, sizeof(size_t));

  std::vector<vid_t> seq(size);
  stream.read((char*)seq.data(), seq.size() * sizeof(vid_t));
  return seq;
}

void writeTextSequence(std::vector<vid_t> const &seq, char const *const filename) {
  std::ofstream stream(filename, std::ios::trunc);
  for (vid_t X : seq)
    stream << X << std::endl;
}

std::vector<vid_t> readTextSequence(char const *const filename) {
  std::vector<vid_t> seq;

  std::ifstream stream(filename);
  vid_t X;
  while(stream >> X)
    seq.push_back(X);

  return seq;
}

void writeSequence(std::vector<vid_t> const &seq, char const *const filename) {
#ifdef USE_BIN_SEQUENCE
  writeBinarySequence(seq, filename);
#else
  writeTextSequence(seq, filename);
#endif
}

std::vector<vid_t> readSequence(char const *const filename) {
#ifdef USE_BIN_SEQUENCE
  return readBinarySequence(filename);
#else
  return readTextSequence(filename);
#endif
}

