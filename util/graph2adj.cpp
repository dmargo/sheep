/*
 * Copyright (c) 2015
 *      The President and Fellows of Harvard College.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

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
    /*
    size_t post_degree = 0;
    for (auto eitr = graph.getEdgeItr(*nitr); !eitr.isEnd(); ++eitr) {
      if (index.at(*eitr) > index.at(*nitr))
        ++post_degree;
    }
    adj << post_degree;
    */
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
