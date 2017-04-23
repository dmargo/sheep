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

#include <algorithm>
#include <chrono>
#include <fstream>
#include <vector>

#include <cmath>
#include <cstdint>

#include <unistd.h>

#include <defs.h>
#ifndef DDUP_GRAPH
  #define DDUP_GRAPH
#endif
#include <graph_wrapper.h>
#include <jnode.h>
#include <jtree.h>
#include <sequence.h>



uint64_t unary_length(uint64_t x) {
  return x; // x - 1 zeros + 1 one
}

uint64_t binary_length(uint64_t x) {
  return std::log2(x) + 1;
}

uint64_t gamma_length(uint64_t x) {
  return unary_length(binary_length(x)) + binary_length(x) - 1;
}

uint64_t sigma_length(uint64_t k, uint64_t x) {
  uint64_t hk = std::log2(x);
  return (hk / k + 1) * (k + 1) - (hk % k == 0 ? 1 : 0);
}

uint64_t delta_length(uint64_t x) {
  return gamma_length(binary_length(x)) + binary_length(x) - 1;
}

uint64_t omega_length(uint64_t x) {
  return x != 1 ? omega_length(binary_length(x) - 1) + binary_length(x) : 1;
}

uint64_t nibble_length(uint64_t x) {
  return ((binary_length(x) / 3) + 1) * 4;
}

void evaluate_gaps(std::vector<size_t> const &gap_count, bool print_cdf = true) {
  uint64_t gamma_sum = 0, gamma_cdf = 0;
  uint64_t sigma_sum = 0, sigma_cdf = 0;
  uint64_t delta_sum = 0, delta_cdf = 0;
  uint64_t omega_sum = 0, omega_cdf = 0;
  uint64_t nibble_sum = 0, nibble_cdf = 0;
  uint64_t u32bit_sum = 0, u32bit_cdf = 0;
  uint64_t u64bit_sum = 0, u64bit_cdf = 0;

  for (size_t x = 1; x < gap_count.size(); ++x) {
    if (gap_count.at(x) != 0) {
      gamma_sum += gap_count.at(x) * gamma_length(x);
      sigma_sum += gap_count.at(x) * sigma_length(2,x);
      delta_sum += gap_count.at(x) * delta_length(x);
      omega_sum += gap_count.at(x) * omega_length(x);
      nibble_sum += gap_count.at(x) * nibble_length(x);
      u32bit_sum += gap_count.at(x) * 32;
      u64bit_sum += gap_count.at(x) * 64;
    }
  }

  if (print_cdf) {
    printf("gap\tgamma\tsigma\tdelta\tomega\tnibble\tu32\tu64\n");
    for (size_t x = 1; x < gap_count.size(); ++x) {
      if (gap_count.at(x) != 0) {
        gamma_cdf += gap_count.at(x) * gamma_length(x);
        sigma_cdf += gap_count.at(x) * sigma_length(2,x);
        delta_cdf += gap_count.at(x) * delta_length(x);
        omega_cdf += gap_count.at(x) * omega_length(x);
        nibble_cdf += gap_count.at(x) * nibble_length(x);
        u32bit_cdf += gap_count.at(x) * 32;
        u64bit_cdf += gap_count.at(x) * 64;

        printf("%zu\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
            x,
            ((double)gamma_cdf)/gamma_sum,
            ((double)sigma_cdf)/sigma_sum,
            ((double)delta_cdf)/delta_sum,
            ((double)omega_cdf)/omega_sum,
            ((double)nibble_cdf)/nibble_sum,
            ((double)u32bit_cdf)/u32bit_sum,
            ((double)u64bit_cdf)/u64bit_sum
        );
      }
    }
  }

  printf("Gamma sum:\t%lu\n", gamma_sum);
  printf("Sigma sum:\t%lu\n", sigma_sum);
  printf("Delta sum:\t%lu\n", delta_sum);
  printf("Omega sum:\t%lu\n", omega_sum);
  printf("Nibble sum:\t%lu\n", nibble_sum);
  printf("32-bit sum:\t%lu\n", u32bit_sum);
  printf("64-bit sum:\t%lu\n", u64bit_sum);
}



enum class pin_mode { BEGIN, NUDGE, MIDDLE, END };

/* XXX sequence_gaps and tree_gaps differ only by a position function and a distance function.
 * It may be wise to eventually template them into one function. */
void sequence_gaps(GraphWrapper const &graph, std::vector<vid_t> const &seq,
    bool forward_only = false, pin_mode mode = pin_mode::BEGIN)
{
  std::vector<size_t> gap_count(seq.size(), 0);

  std::vector<size_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1);
  for (size_t i = 0; i != seq.size(); ++i)
    pos.at(seq.at(i)) = i;
  auto gap_between = [](size_t const lesser_pos, size_t const greater_pos) -> size_t
  {
    assert(lesser_pos <= greater_pos);
    return greater_pos - lesser_pos;
  };

  std::vector<size_t> buf;
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
    vid_t const X = *nitr;
    size_t const X_pos = pos.at(X);

    if (buf.size() < graph.getDeg(X)) buf.resize(graph.getDeg(X));
    auto buf_itr = buf.begin();

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      size_t const Y_pos = pos.at(Y);

      if (X_pos != Y_pos || (forward_only && X_pos < Y_pos))
        *(buf_itr++) = Y_pos;
    }
    if (buf_itr != buf.begin()) {
      std::sort(buf.begin(), buf_itr);
      size_t const buf_size = std::distance(buf.begin(), buf_itr);
      size_t const pinned_pos = 
        mode == pin_mode::BEGIN ? buf.at(0) :
        mode == pin_mode::NUDGE ? buf.at(buf_size == 1 ? 0 : 1) :
        mode == pin_mode::MIDDLE ? buf.at(buf_size/2) :
        buf.at(buf_size-1); // mode == pin_mode::END
      size_t const initial_gap = pinned_pos <= X_pos ?
        gap_between(pinned_pos, X_pos) :
        gap_between(X_pos, pinned_pos);
      assert(initial_gap != 0);
      gap_count.at(initial_gap) += 1;
      for (size_t i = 1; i != buf_size; ++i) {
        size_t const this_gap = gap_between(buf.at(i-1), buf.at(i));
        assert(this_gap != 0);
        gap_count.at(this_gap) += 1;
      }
    }
  }

  assert(gap_count.at(0) == 0);
  evaluate_gaps(gap_count, false);
}

void tree_gaps(GraphWrapper const &graph, JTree const &tree,
    bool forward_only = false, pin_mode mode = pin_mode::BEGIN)
{
  std::vector<size_t> gap_count(tree.jnodes.size());

  std::vector<size_t> depth(tree.jnodes.size(), 0);
  for (jnid_t i = tree.jnodes.size() - 1; i != (jnid_t)-1; --i) {
    jnid_t const parent = tree.jnodes.parent(i);
    depth.at(i) = parent != INVALID_JNID ? depth.at(parent) + 1 : 0;
  }

  std::vector<jnid_t> buf;
  for (auto nitr = graph.getNodeItr(); !nitr.isEnd(); ++nitr) {
    vid_t const X = *nitr;
    jnid_t const X_pos = tree.vid2jnid(X);

    if (buf.size() < graph.getDeg(X)) buf.resize(graph.getDeg(X));
    auto buf_itr = buf.begin();

    for (auto eitr = graph.getEdgeItr(X); !eitr.isEnd(); ++eitr) {
      vid_t const Y = *eitr;
      jnid_t const Y_pos = tree.vid2jnid(Y);

      if (X_pos != Y_pos || (forward_only && X_pos < Y_pos))
        *(buf_itr++) = Y_pos;
    }
    if (buf_itr != buf.begin()) {
      auto gap_between = [&depth,X_pos](jnid_t const lesser_pos, jnid_t const greater_pos) -> size_t
      {
        assert(lesser_pos <= greater_pos);
        if (greater_pos <= X_pos)       // Both are less
          return greater_pos - lesser_pos;
        else if (X_pos <= lesser_pos)   // Both are greater
          return depth.at(lesser_pos) - depth.at(greater_pos);
        else {                          assert(lesser_pos < X_pos && X_pos < greater_pos);
          return (X_pos - lesser_pos) + (depth.at(X_pos) - depth.at(greater_pos));
        }
      };

      std::sort(buf.begin(), buf_itr);
      size_t const buf_size = std::distance(buf.begin(), buf_itr);
      size_t const pinned_pos = 
        mode == pin_mode::BEGIN ? buf.at(0) :
        mode == pin_mode::NUDGE ? buf.at(buf_size == 1 ? 0 : 1) :
        mode == pin_mode::MIDDLE ? buf.at(buf_size/2) :
        buf.at(buf_size-1); // mode == pin_mode::END
      size_t const initial_gap = pinned_pos <= X_pos ?
        gap_between(pinned_pos, X_pos) :
        gap_between(X_pos, pinned_pos);
      assert(initial_gap != 0);
      gap_count.at(initial_gap) += 1;
      for (size_t i = 1; i != buf_size; ++i) {
        size_t const this_gap = gap_between(buf.at(i-1), buf.at(i));
        assert(this_gap != 0);
        gap_count.at(this_gap) += 1;
      }
    }
  }

  assert(gap_count.at(0) == 0);
  evaluate_gaps(gap_count, false);
}

void tree_cost(JTree const &tree) {
  std::vector<size_t> gap_count(tree.jnodes.size());
  for (jnid_t id = 0; id != tree.jnodes.size(); ++id) {
    if (tree.jnodes.parent(id) == INVALID_JNID)
      gap_count.at(1) += 1;
    else
      gap_count.at(tree.jnodes.parent(id) - id) += 1;
  }

  evaluate_gaps(gap_count, false);
}



int main(int argc, char* argv[]) {
  bool undirected = true;
  bool forward_only = false;
  pin_mode mode = pin_mode::BEGIN;

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "dfnme")) != -1) {
    switch (opt) {
      case 'd':
        undirected = false;
        break;
      case 'f':
        forward_only = true;
        break;
      case 'n':
        mode = pin_mode::NUDGE;
        break;
      case 'm':
        mode = pin_mode::MIDDLE;
        break;
      case 'e':
        mode = pin_mode::END;
        break;
      case '?':
        printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if (optind >= argc) {
    printf("USAGE: gaps graph [seq] [tree]\n");
    return 1;
  }
  char const *const graph_filename = argv[optind];
  char const *const seq_filename = optind + 1 < argc ? argv[optind + 1] : "";
  char const *const tree_filename = optind + 2 < argc ? argv[optind + 2] : "";
   
  auto start_point = std::chrono::steady_clock::now();

  GraphWrapper graph(graph_filename, 0,0, undirected);
  printf("Nodes:%zu Edges:%zu\n", graph.getNodes(), graph.getEdges());

  auto load_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  //printf("Loaded in: %lums\n\n", load_duration.count());

  std::vector<vid_t> seq = strcmp(seq_filename, "") != 0 ?
    readSequence(seq_filename) :
    identitySequence(graph);

  if (strcmp(tree_filename, "") == 0) {
    sequence_gaps(graph, seq, forward_only, mode);
  }
  else {
    JTree tree(seq, tree_filename);
    tree_gaps(graph, tree, forward_only, mode);
    tree_cost(tree);
  }
  
  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      (std::chrono::steady_clock::now() - start_point) - load_duration);
  //printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
