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
#include <jnode.h>

int main(int argc, char* argv[]) {

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

  size_t edge_count = 0;
  std::vector<size_t> edge_width(jnodes.size(), 0);
  std::vector<size_t> subt(jnodes.size(), 1);
  std::vector<size_t> supr(jnodes.size(), 1);

  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    edge_width.at(id) += jnodes.pst_weight(id);
    if (jnodes.parent(id) != INVALID_JNID) {
      ++edge_count;
      edge_width.at(jnodes.parent(id)) += edge_width.at(id) - jnodes.pre_weight(id);
      subt.at(jnodes.parent(id)) += subt.at(id);
    }
  }
  for (jnid_t id = jnodes.size() - 1; id != (jnid_t) -1; --id) {
    if (jnodes.parent(id) != INVALID_JNID) {
      supr.at(id) += supr.at(jnodes.parent(id));
    }
  }

  std::ofstream adj(adj_filename);
  adj << jnodes.size() << ' ' << edge_count << " 011" << std::endl;
  for (jnid_t id = 0; id != jnodes.size(); ++id) {
    adj << 1; // vertex weight

    jnid_t const par_id = jnodes.parent(id);
    if (par_id != INVALID_JNID)
      adj << ' ' << par_id + 1 << ' ' <<
        (std::min(subt.at(id),edge_width.at(id)) + std::min(supr.at(par_id),edge_width.at(id)));

    for (auto kid_id : jnodes.kids(id))
      adj << ' ' << kid_id + 1 << ' ' <<
        (std::min(subt.at(kid_id),edge_width.at(kid_id)) + std::min(supr.at(id),edge_width.at(kid_id)));

    adj << std::endl;
  }
  
  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
