
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <vector>

#include <defs.h>
#include <jnode.h>

int main(int argc, char* argv[]) {

  if (optind + 1 >= argc) {
    printf("USAGE: graph2dot input_graph output_graph\n");
    return 1;
  }
  char const *const tree_filename = argv[optind];
  char const *const dot_filename = argv[optind + 1];
   
  auto start_point = std::chrono::steady_clock::now();

  JNodeTable jnodes(tree_filename);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Loaded in: %lums\n\n", load_duration.count());

  std::ofstream dot(dot_filename);
  dot << "digraph {" << std::endl;
  for (jnid_t id = jnodes.size() - 1; id != (jnid_t)-1; --id) {
    dot << "\t" << id;
    if (jnodes.parent(id) != INVALID_JNID)
      dot << " -> " << jnodes.parent(id);
    dot << std::endl;
  }
  dot << "}" << std::endl;
  
  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
