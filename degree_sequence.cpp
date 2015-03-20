#include <chrono>
#include <unistd.h>

#include <sequence.h>
#include <stdafx.h>

int main(int argc, char* argv[]) {

  char const *output_filename = "out.seq";

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "o:")) != -1) {
    switch (opt) {
      case 'o':
        output_filename = optarg;
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

  if (optind >= argc) {
    printf("USAGE: degree_sequence [options] graph_file");
    return 1;
  }
  
  auto start_point = std::chrono::steady_clock::now();

  std::vector<vid_t> seq = fileSequence(argv[optind]);
  writeSequence(seq, output_filename);

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point));
  printf("Sorted in: %lums\n", run_duration.count());

  return 0;
}
