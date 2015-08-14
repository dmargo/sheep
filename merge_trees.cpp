
#include <chrono>
#include <unistd.h>
#include <vector>

#include <defs.h>
#include <jnode.h>

int main(int argc, char* argv[]) {
  char const *output_filename = "";

  bool verbose = false;
  bool make_kids = false;
  bool do_faqs = false;

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "o:vkf")) != -1) {
    switch (opt) {
      case 'o':
        output_filename = optarg;
        break;
      case 'v':
        verbose = !verbose;
        break;
      case 'k':
        make_kids = !make_kids;
        break;
      case 'f':
        do_faqs = !do_faqs;
        break;
      case '?':
        if (optopt == 'o')
          printf("Option -%c requires a string.\n", optopt);
        else
          printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if (optind + 1 >= argc) {
    printf("USAGE: merge_trees [options ...] first.tree second.tree\n");
    return 1;
  }
   
  auto start_point = std::chrono::steady_clock::now();

  JNodeTable lhs(argv[optind]);
  JNodeTable rhs(argv[optind + 1]);

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  if (verbose) printf("Loaded in: %lums\n", load_duration.count());

  JNodeTable jnodes = strcmp(output_filename, "") == 0 ?
    JNodeTable(lhs.size(), make_kids, 0) :
    JNodeTable(output_filename, lhs.size(), make_kids, 0);
  jnodes.merge(lhs, rhs, make_kids);

  auto build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  if (verbose) printf("Built in: %lums\n", build_duration.count());

  if (do_faqs) {
    JNodeTable::Facts faq = jnodes.getFacts();
    faq.print();
  }

  return 0;
}
