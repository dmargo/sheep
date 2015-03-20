#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>

#include <mpi.h>
#include <unistd.h>

#include <graph_wrapper.h>
#include <jtree.h>
#include <sequence.h>
#include <stdafx.h>

int main(int argc, char* argv[]) {
  bool use_mpi_sort = false;
  bool use_mpi_reduce = false;

  size_t part = 0;
  size_t num_parts = 0;
  char const *sequence_filename = "";
  char const *output_filename = "";

  auto jopts = JTree::Options();
  bool do_faqs = false;
  bool do_print = false;
  bool do_validate = false;

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "irl:s:o:vkejm:w:xfdpc")) != -1) {
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
      case 'p':
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
  
    // If using MPI and mmap together, make sure maps have different names.
    if (!use_mpi_reduce && strcmp(output_filename, "") != 0) {
      char *const tmp_filename = (char*)malloc(strlen(output_filename) + 9);
      sprintf(tmp_filename, "%s%02dr0.tre", output_filename, rank);
      output_filename = tmp_filename;
    }
  }
  bool const is_leader = ((use_mpi_sort || use_mpi_reduce) && part == 1) ||
                         (!(use_mpi_sort || use_mpi_reduce) && strcmp(sequence_filename, "") == 0);

  if (jopts.verbose) printf("Loading %s...\n", graph_filename);
  GraphWrapper graph(graph_filename, part, num_parts);
  if (jopts.verbose) printf("Nodes:%u Edges:%zu\n", graph.getNodes(), graph.getEdges());

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  if (is_leader) printf("Loaded graph in: %f seconds\n", load_duration.count() / 1000.0);
  
  std::vector<vid_t> seq = use_mpi_sort ?
    mpiSequence(graph) : strcmp(sequence_filename, "") != 0 ?
    readSequence(sequence_filename) :
    degreeSequence(graph);

  auto sort_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point) - load_duration;
  if (is_leader && (use_mpi_sort || strcmp(sequence_filename, "") == 0))
    printf("Sorted in: %f seconds\n", sort_duration.count() / 1000.0);

  JTree tree = !use_mpi_reduce && (strcmp(output_filename, "") != 0) ?
    JTree(graph, seq, output_filename, jopts) :
    JTree(graph, seq, jopts);

  auto map_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point) - sort_duration - load_duration;
  if (is_leader) printf("Mapped in: %f seconds\n", map_duration.count() / 1000.0);

  if (use_mpi_reduce) {
    tree.jnodes.mpi_merge();

    auto reduce_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_point) - map_duration - sort_duration - load_duration;
    if (is_leader) printf("Reduced in: %f seconds\n", reduce_duration.count() / 1000.0);
  }

  if (use_mpi_sort || use_mpi_reduce) {
    MPI_Finalize();
    if (use_mpi_sort && part == 1 && strcmp(sequence_filename, "") != 0)
      writeSequence(seq, sequence_filename);
    if (use_mpi_reduce && part == 1 && strcmp(output_filename, "") != 0)
      tree.jnodes.save(output_filename);
  }

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
