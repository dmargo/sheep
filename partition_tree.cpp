
#include <chrono>
#include <unistd.h>
#include <vector>

#include <defs.h>
#include <graph_wrapper.h>
#include <jnode.h>
#include <partition.h>
#include <sequence.h>

int main(int argc, char* argv[]) {

  bool verbose = true;
  bool do_faqs = false;

  double balance_factor = 1.03;
  bool vtx_weight = false;
  bool pst_weight = false;
  bool pre_weight = false;

  char const *graph_filename = "";
  char const *output_filename = "";

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "vfb:xdug:o:")) != -1) {
    switch (opt) {
      case 'v':
        verbose = !verbose;
        break;
      case 'f':
        do_faqs = !do_faqs;
        break;
      case 'b':
        balance_factor = atof(optarg);
        break;
      case 'x':
        vtx_weight = true;
        break;
      case 'd':
        pst_weight = true;
        break;
      case 'u':
        pre_weight = true;
        break;
      case 'g':
        graph_filename = optarg;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case '?':
        if (optopt == 'k')
          printf("Option -%c requires a long long.\n", optopt);
        else if (optopt == 'b')
          printf("Option -%c requires a double.\n", optopt);
        else if (optopt == 'g' || optopt == 'o')
          printf("Option -%c requires a string.\n", optopt);
        else
          printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }
  if (!(vtx_weight || pst_weight || pre_weight))
    pst_weight = true;

  if (optind + 2 >= argc) {
    printf("USAGE: partition_tree [options] input_sequence input_tree parts [parts...]\n");
    return 1;
  }
   
  auto start_point = std::chrono::steady_clock::now();

  JNodeTable jnodes(argv[optind + 1]);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  if (verbose) printf("Loaded tree in: %f seconds\n", load_duration.count() / 1000.0);

  if (do_faqs)
    jnodes.getFacts().print();

  /* SIMPLE PARTITIONING */
  if (strcmp(graph_filename, "") == 0) {
    std::vector<vid_t> seq = readSequence(argv[optind]);
    for (int i = optind + 2; i != argc; ++i) {
      short const num_parts = atoi(argv[optind + 2]);
      Partition part(seq, jnodes, num_parts, balance_factor, vtx_weight, pst_weight, pre_weight);
      part.print();
    }
  }
  /* PARTITIONING AND EVALUATION */
  else if (strcmp(output_filename, "") == 0) {
    GraphWrapper graph(graph_filename);
    std::vector<vid_t> seq = strcmp(argv[optind], "-") == 0 ?
      degreeSequence(graph) :
      readSequence(argv[optind]);

    for (int i = optind + 2; i != argc; ++i) {
      short const num_parts = atoi(argv[i]);

      auto partition_start = std::chrono::steady_clock::now();

      Partition part(seq, jnodes, num_parts, balance_factor, vtx_weight, pst_weight, pre_weight);

      auto partition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - partition_start);
      if (verbose) printf("Partitioning took: %f seconds\n", partition_duration.count() / 1000.0);

      part.print();
      part.evaluate(graph, seq);
    }
  }
  /* PARTITIONING AND I/O */
  else {
    std::vector<vid_t> seq = strcmp(argv[optind], "-") == 0 ?
      fileSequence(graph_filename) :
      readSequence(argv[optind]);

    short const num_parts = atoi(argv[optind + 2]);

    auto partition_start = std::chrono::steady_clock::now();

    Partition part(seq, jnodes, num_parts, balance_factor, vtx_weight, pst_weight, pre_weight);

    auto partition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::steady_clock::now() - partition_start);
    if (verbose) printf("Partitioning took: %f seconds\n", partition_duration.count() / 1000.0);

    part.print();
    part.writePartitionedGraph(graph_filename, seq, output_filename);
  }

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  if (verbose) printf("Finished in: %f seconds\n", run_duration.count() / 1000.0);

  return 0;
}
