#pragma once

#include <algorithm>
#include <fstream>
#include <vector>

#include <mpi.h>
#include <parallel/algorithm>

#include "defs.h"

/* SEQEUENCE CONSTRUCTORS */
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
  vid_t max_vid = 0;
  vid_t local_max = graph.getMaxVid();
  MPI_Allreduce((void*)&local_max, (void*)&max_vid, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

  vid_t *degree = (vid_t*)malloc((max_vid + 1) * sizeof(vid_t));
  vid_t *local_degree = (vid_t*)calloc(max_vid + 1, sizeof(vid_t));
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr)
    local_degree[*nitr] = graph.getDeg(*nitr);
  MPI_Allreduce((void*)local_degree, (void*)degree, max_vid, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
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
  std::vector<vid_t> degree;

  if (strcmp(".net", filename + strlen(filename) - 4) == 0) {
    vid_t X;
    std::ifstream stream(filename);
    while(stream >> X) {
      if (degree.size() < X + 1) degree.resize(X + 1);
      degree[X] += 1;
    }
  } else if (strcmp(".dat", filename + strlen(filename) - 4) == 0) {
    struct xs1 {
	    unsigned tail;
	    unsigned head;
	    float weight;
    };
    xs1 buf;
    std::ifstream stream(filename, std::ios::binary);
    while (!stream.eof()) {
      stream.read((char*)&buf, sizeof(xs1));
      size_t const required_size = std::max(buf.tail, buf.head) + 1;
      if (degree.size() < required_size) degree.resize(required_size);
      degree[buf.tail] += 1;
      degree[buf.head] += 1;
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

