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
#include <vector>

#include <unistd.h>

#include <defs.h>
#include <jnode.h>
#include <jtree.h>
#include <readerwriter.h>
#include <sequence.h>

/* This duplicates code in lib/partition.cpp;
 * that code should probably be factored out into a generic function. */
template <typename ReaderType = SNAPReader, typename WriterType = SNAPWriter>
void write_isomorphic_graph(
    char const *const input_filename, std::vector<vid_t> const &seq,
    char const *const output_filename, bool forward_only = false)
{
  /* This code is also extremely generic and probably belongs in sequence.h */
  std::vector<jnid_t> pos(*std::max_element(seq.cbegin(), seq.cend()) + 1, INVALID_JNID);
  for (jnid_t i = 0; i != seq.size(); ++i)
    pos[seq[i]] = i;

  vid_t X,Y;
  ReaderType reader(input_filename);
  WriterType writer(output_filename);
  while(reader.read(X,Y)) {
    jnid_t X_pos = pos.at(X);
    jnid_t Y_pos = pos.at(Y);
    if (!forward_only || X_pos < Y_pos)
      writer.write(X_pos,Y_pos);
  }
}

template <typename ReaderType = SNAPReader, typename WriterType = SNAPWriter>
void write_transformed_graph(
    char const *const input_filename, std::vector<vid_t> const &seq, JTree const &tree,
    char const *const output_filename)
{
  std::vector<size_t> depth(tree.jnodes.size(), 0);
  for (jnid_t i = tree.jnodes.size() - 1; i != (jnid_t)-1; --i) {
    jnid_t const parent = tree.jnodes.parent(i);
    depth.at(i) = parent != INVALID_JNID ? depth.at(parent) + 1 : 0;
  }

  vid_t X,Y;
  ReaderType reader(input_filename);
  WriterType writer(output_filename);
  while(reader.read(X,Y)) {
    jnid_t X_pos = tree.vid2jnid(X);
    jnid_t Y_pos = tree.vid2jnid(Y);
    if (X_pos < Y_pos) {
      size_t gap = depth.at(X_pos) - depth.at(Y_pos);
      writer.write(X_pos,X_pos + gap);
    }
    else {
      writer.write(X_pos,Y_pos);
    }
  }
}

int main(int argc, char* argv[]) {
  bool forward_only = false;
  char const *output_filename = "out"; // Should really start using std::cout

  opterr = 0;
  int opt;
  while ((opt = getopt(argc, argv, "fo:")) != -1) {
    switch (opt) {
      case 'f':
        forward_only = true;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case '?':
        printf("Unknown option character '\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if (optind + 1 >= argc) {
    printf("USAGE: transform graph seq [tree]\n");
    return 1;
  }
  char const *const graph_filename = argv[optind];
  char const *const seq_filename = argv[optind + 1];
  char const *const tree_filename = optind + 2 < argc ? argv[optind + 2] : "";
   
  auto start_point = std::chrono::steady_clock::now();

  std::vector<vid_t> seq = readSequence(seq_filename);
  if (strcmp(tree_filename, "") == 0) {
    write_isomorphic_graph(graph_filename, seq, output_filename, forward_only);
  }
  else {
    JTree tree(seq, tree_filename);
    write_transformed_graph(graph_filename, seq, tree, output_filename);
  }

  auto run_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_point);
  //printf("Finished in: %lums\n", run_duration.count());

  return 0;
}
