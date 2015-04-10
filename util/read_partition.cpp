
#include <chrono>
#include <unistd.h>
#include <vector>

#include <graph_wrapper.h>
#include <partition.h>
#include <sequence.h>
#include <stdafx.h>

int main(int argc, char* argv[]) {

  if (argc < 3) {
    printf("USAGE: read_partition graph partition [partition...]\n");
    return 1;
  }
   
  auto start_point = std::chrono::steady_clock::now();

  GraphWrapper graph(argv[1]);
  std::vector<vid_t> seq = degreeSequence(graph);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Loaded in: %lums\n", load_duration.count());
  printf("Nodes:%u Edges:%zu\n", graph.getNodes(), graph.getEdges());

  for (int i = 2; i < argc; ++i) {
    Partition part(seq, argv[i]);
    part.evaluate(graph);
  }

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
