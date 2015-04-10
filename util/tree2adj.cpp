
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <vector>

#include <jnode.h>
#include <stdafx.h>

int main(int argc, char* argv[]) {

#if 0
  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "fk:g:")) != -1) {
    switch (opt) {
      case 'f':
        do_faqs = !do_faqs;
        break;
      case 'k':
        num_parts = atoll(optarg);
        break;
      case 'g':
        graph_filename = optarg;
        break;
      case '?':
        if (optopt == 'c')
          printf("Option -%c requires a long long.\n", optopt);
        else if (optopt == 'g')
          printf("Option -%c requires a string.\n", optopt);
        else
          printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }
#endif

  if (optind + 1 >= argc) {
    printf("USAGE: tree2adj input_tree output_graph\n");
    return 1;
  }
  char const *const tree_filename = argv[optind];
  char const *const adj_filename = argv[optind + 1];
   
  auto start_point = std::chrono::steady_clock::now();

  JNodeTable jnodes(tree_filename);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  printf("Loaded in: %lums\n\n", load_duration.count());

  size_t edge_cnt = 0;
  std::vector<size_t> subt(jnodes.size(), 1);
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    if (jnodes.parent(id) != INVALID_JNID) {
      ++edge_cnt;
      subt.at(jnodes.parent(id)) += subt.at(id);
    }
  }

  std::vector<size_t> supr(jnodes.size(), 1);
  for (jnid_t id = jnodes.size() - 1; id != (jnid_t) -1; --id) {
    if (jnodes.parent(id) != INVALID_JNID) {
      supr.at(id) += supr.at(jnodes.parent(id));
    }
  }

  std::ofstream adj(adj_filename);
  adj << jnodes.size() << ' ' << edge_cnt << " 011" << std::endl;
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    adj << 1; // vertex weight

    jnid_t const par_id = jnodes.parent(id);
    if (par_id != INVALID_JNID)
      adj << ' ' << par_id + 1 << ' ' << (subt.at(id) + supr.at(par_id));
     
    for (auto itr = jnodes.kids(id).cbegin(); itr != jnodes.kids(id).cend(); ++itr) {
      jnid_t const kid_id = *itr;
      adj << ' ' << kid_id + 1 << ' ' << (subt.at(kid_id) + supr.at(id));
    }

    adj << std::endl;
  }
  
  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
