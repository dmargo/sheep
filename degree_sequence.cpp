
#include <chrono>

#include <sequence.h>
#include <stdafx.h>

int main(int argc, char* argv[]) {
  if (argc != 3) {
    printf("USAGE: degree_sequence graph_file output_file");
    return 1;
  }
  
  auto start_point = std::chrono::steady_clock::now();

  std::vector<vid_t> seq = fileSequence(argv[1]);
  writeSequence(seq, argv[2]);

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point));
  printf("Sorted in: %lums\n", run_duration.count());

  return 0;
}
