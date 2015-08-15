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

#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>

#include <mpi.h>
#include <unistd.h>

#include <defs.h>
#include <graph_wrapper.h>
#include <jtree.h>
#include <partition.h>
#include <sequence.h>

int main(int argc, char* argv[]) {
  bool use_mpi_sort = false;
  bool use_mpi_reduce = false;

  size_t part = 0;
  size_t num_parts = 0;
  size_t partitions = 0;
  char const *sequence_filename = "";
  char const *output_filename = "";

  auto jopts = JTree::Options();
  bool do_faqs = false;
  bool do_print = false;
  bool do_validate = false;

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "irl:p:s:o:vkejm:w:xfdtc")) != -1) {
    switch (opt) {
      case 'i':
        use_mpi_sort = !use_mpi_sort;
        break;
      case 'r':
        use_mpi_reduce = !use_mpi_reduce;
        break;
      case 'l':
        part=atoll(strtok(optarg, "/"));
        num_parts=atoll(strtok(nullptr, "/"));
        assert(strtok(nullptr, "/") == nullptr);
        break;
      case 'p':
        partitions=atoll(optarg);
        break;
      case 's':
        sequence_filename = optarg;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'v':
        jopts.verbose = !jopts.verbose;
        break;
      case 'k':
        jopts.make_kids = !jopts.make_kids;
        break;
      case 'e':
        jopts.make_pst = !jopts.make_pst;
        break;
      case 'j':
        jopts.make_jxn = !jopts.make_jxn;
        break;
      case 'm':
        jopts.memory_limit = atoll(optarg) * MEGA;
        break;
      case 'w':
        jopts.width_limit = atoll(optarg);
        break;
      case 'x':
        jopts.find_max_width = !jopts.find_max_width;
        break;
      case 'f':
        do_faqs = !do_faqs;
        break;
      case 't':
        do_print = !do_print;
        break;
      case 'c':
        do_validate = !do_validate;
        break;
      case '?':
        if (optopt == 's' || optopt == 'o')
          printf("Option -%c requires a string.\n", optopt);
        else if (optopt == 'm' || optopt == 'w')
          printf("Option -%c requires a long long.\n", optopt);
        else
          printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if (optind >= argc) {
    printf("USAGE: graph2tree input_graph [options ...]\n");
    return 1;
  }
  char const *const graph_filename = argv[optind];
  
  auto start_point = std::chrono::steady_clock::now();

  if (use_mpi_sort || use_mpi_reduce) {
    MPI_Init(nullptr, nullptr);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    part = rank + 1;

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    num_parts = size;
  
    // If using MPI and mmap, make sure maps have different names.
    if (!use_mpi_reduce && strcmp(output_filename, "") != 0) {
      char *const tmp_filename = (char*)malloc(strlen(output_filename) + 9);
      sprintf(tmp_filename, "%s%02dr0.tre", output_filename, rank);
      output_filename = tmp_filename;
    }
    // If using MPI to output partitions, make sure outputs have different names.
    else if (use_mpi_reduce && partitions != 0 && strcmp(output_filename, "") != 0) {
      char *const tmp_filename = (char*)malloc(strlen(output_filename) + 9);
      sprintf(tmp_filename, "%s-w%04d-p", output_filename, rank);
      output_filename = tmp_filename;
    }
  }
  bool const is_leader = ((use_mpi_sort || use_mpi_reduce) && part == 1) ||
                         (!(use_mpi_sort || use_mpi_reduce) && strcmp(sequence_filename, "") == 0);

  if (jopts.verbose) printf("Loading %s...\n", graph_filename);
  GraphWrapper graph(graph_filename, part, num_parts);
  if (jopts.verbose) printf("Nodes:%zu Edges:%zu\n", graph.getNodes(), graph.getEdges());

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  if (is_leader) printf("Loaded graph in: %f seconds\n", load_duration.count() / 1000.0);
  
  std::vector<vid_t> seq =
    use_mpi_sort ?
      mpiSequence(graph) :
    strcmp(sequence_filename, "") != 0 ?
      readSequence(sequence_filename) :
    //else
      degreeSequence(graph);

  if (use_mpi_sort && part == 1 && strcmp(sequence_filename, "") != 0)
    writeSequence(seq, sequence_filename);

  auto sort_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point) - load_duration;
  if (is_leader && (use_mpi_sort || strcmp(sequence_filename, "") == 0))
    printf("Sorted in: %f seconds\n", sort_duration.count() / 1000.0);

  JTree tree =
    !use_mpi_reduce && strcmp(output_filename, "") != 0 && partitions == 0 ?
      JTree(graph, seq, output_filename, jopts) :
    //else
      JTree(graph, seq, jopts);

  auto map_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point) - sort_duration - load_duration;
  if (is_leader) printf("Mapped in: %f seconds\n", map_duration.count() / 1000.0);

  if (use_mpi_reduce) {
    tree.jnodes.mpi_merge(jopts.make_kids);

    auto reduce_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_point) - map_duration - sort_duration - load_duration;
    if (is_leader) printf("Reduced in: %f seconds\n", reduce_duration.count() / 1000.0);
  }

  if (partitions != 0) {
    //XXX Unfortunate hack; Partition requires kids, which are normally made on-load.
    if (!use_mpi_reduce || part == 1)
      tree.jnodes.makeKids();
    Partition p = !use_mpi_reduce || part == 1 ?
      Partition(seq, tree.jnodes, partitions) : Partition();
    if (use_mpi_reduce)
      p.mpi_sync();

    if (strcmp(output_filename, "") != 0)
      p.writePartitionedGraph(graph, seq, output_filename);
    else if (is_leader)
      p.print();
  }
  else if (use_mpi_reduce && part == 1 && strcmp(output_filename, "") != 0)
    tree.jnodes.save(output_filename);

  if (use_mpi_sort || use_mpi_reduce)
    MPI_Finalize();

  auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point));
  if (jopts.verbose) printf("Built in: %f seconds\n", build_duration.count() / 1000.0);

  if (do_faqs)
    tree.jnodes.getFacts().print();
  if (do_print)
    tree.print();
  if (do_validate) {
    if (tree.isValid(graph, seq, jopts))
      printf("Tree is valid.\n");
    else
      printf("ERROR: Tree is not valid.\n");
  }

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point));
  if (jopts.verbose) printf("Finished in: %f seconds\n", run_duration.count() / 1000.0);

  return 0;
}

