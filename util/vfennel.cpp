
#include <chrono>
#include <unistd.h>
#include <vector>

#include <defs.h>
#include <graph_wrapper.h>
#include <partition.h>
#include <sequence.h>

int main(int argc, char* argv[]) {

  if (argc < 3) {
    printf("USAGE: vfennel graph parts [parts...]\n");
    return 1;
  }
   
  auto start_point = std::chrono::steady_clock::now();

  GraphWrapper graph(argv[1]);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Loaded in: %lums\n", load_duration.count());
  printf("Nodes:%u Edges:%zu\n", graph.getNodes(), graph.getEdges());

  std::vector<vid_t> seq = defaultSequence(graph);

  for (int i = 2; i != argc; ++i) {
    short const num_parts = atoi(argv[i]);

    auto partition_start = std::chrono::steady_clock::now();

    Partition part(graph, seq, num_parts);
    part.print();

    auto partition_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::steady_clock::now() - partition_start);
    printf("Partitioning took: %lums\n", partition_duration.count());

    part.evaluate(graph);
  }

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
