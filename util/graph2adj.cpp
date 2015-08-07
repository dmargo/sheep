
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <vector>

#include <defs.h>
#include <graph_wrapper.h>
#include <sequence.h>

int main(int argc, char* argv[]) {
  if (optind + 1 >= argc) {
    printf("USAGE: graph2adj input_graph output_graph\n");
    return 1;
  }
  char const *const graph_filename = argv[optind];
  char const *const adj_filename = argv[optind + 1];
   
  auto start_point = std::chrono::steady_clock::now();

  GraphWrapper graph(graph_filename);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Loaded in: %lums\n\n", load_duration.count());

  std::vector<vid_t> seq = degreeSequence(graph);

  std::vector<vid_t> index(*std::max_element(seq.cbegin(), seq.cend()) + 1);
  for (size_t i = 0; i != seq.size(); ++i)
    index.at(seq.at(i)) = i + 1;

  //Have to compute this because LLAMA counts self-edges.
  size_t edge_cnt = 0;
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr)
    for (auto eitr = graph.getEdgeItr(*nitr); !eitr.isEnd(); ++eitr)
      if (*eitr > *nitr)
        ++edge_cnt;

  std::ofstream adj(adj_filename);
  adj << graph.getNodes() << ' ' << edge_cnt << ' ' << "010" << std::endl;
  for (auto nitr = seq.cbegin(); nitr != seq.cend(); ++nitr) {
    size_t post_degree = 0;
    for (auto eitr = graph.getEdgeItr(*nitr); !eitr.isEnd(); ++eitr) {
      if (index.at(*eitr) > index.at(*nitr))
        ++post_degree;
    }
    adj << graph.getDeg(*nitr);
    for (auto eitr = graph.getEdgeItr(*nitr); !eitr.isEnd(); ++eitr) {
      if (*eitr != *nitr)
        adj << ' ' << index.at(*eitr);
    }
    adj << std::endl;
  }
  
  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
