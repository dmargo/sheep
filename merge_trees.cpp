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

