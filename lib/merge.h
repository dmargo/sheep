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

#pragma once

#include <queue>

#include "defs.h"
#include "jdata.h"

struct SortedRange {
  vid_t const *itr;
  vid_t const *end;

  inline SortedRange(vid_t const *begin, vid_t const *e) : itr(begin), end(e) {}
  inline bool operator>(SortedRange const &rhs) const { return (*itr > *rhs.itr); }
  inline size_t size() const { return std::distance(itr,end); }
};

bool balance_line_merge(JData<vid_t> &new_data, size_t max_len,
    std::vector<SortedRange> &kid_itrs, vid_t Xclude)
{
  for (auto end = kid_itrs.end(); end != kid_itrs.begin();) {
    // Find range of minimums (linear search)
    vid_t min = std::numeric_limits<vid_t>::max();
    std::vector<SortedRange>::iterator frst, last;
    for (auto cur = kid_itrs.begin(); cur != end; cur++) {
      if (*cur->itr < min) {
        min = *cur->itr;
        frst = last = cur;
      } else if (*cur->itr == min) {
        last = cur;
      }
    }

    // Insert min; always need to exclude one elimination element.
    if (min != Xclude) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(min);
    }
    
    // Increment min (linear search)
    for (; frst != last + 1 && frst != end; frst++)
      if (*frst->itr == min)
        if (++(frst->itr) == frst->end)
          std::iter_swap(--end, frst--);
  }
  return true;
}

/* heap merge may outscale balance-line merge for something like kid_itrs.size() > 32.
 * The likelihood of this case increases with graph density.
 * XXX: TBD whether this is significant for Big graphs.
 * XXX: I've seen at least one case where it made a HUGE difference. Needs to be benchmarked.
 * does priority_queue virtualize std::greater calls? */
bool heap_merge(JData<vid_t> &new_data, size_t max_len,
    std::vector<SortedRange> &kid_itrs, vid_t Xclude)
{
  std::priority_queue< SortedRange, std::vector<SortedRange>, std::greater<SortedRange> >
      q(std::greater<SortedRange>(), std::move(kid_itrs));
  vid_t prev = INVALID_VID;

  while (!q.empty()) {
    SortedRange cur = q.top();
    q.pop();

    vid_t const min = *cur.itr++;
    if (min != Xclude && min != prev) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(min);
      prev = min;
    }

    if (cur.itr != cur.end)
      q.push(cur);
  }
  return true;
}

bool asymmetric_merge(JData<vid_t> &new_data, size_t max_len,
    std::vector<SortedRange> &kid_itrs, vid_t Xclude)
{
  assert(kid_itrs.size() == 2);
  SortedRange big = kid_itrs[0];
  SortedRange small = kid_itrs[1];
  if (big.size() < small.size())
    std::swap(big,small);
  if (big.size() < small.size() * 8)
    return balance_line_merge(new_data, max_len, kid_itrs, Xclude);

  for (; small.itr != small.end; small.itr++) {
    auto *big_middle = std::lower_bound(big.itr, big.end, *small.itr);

    for (; big.itr != big_middle; big.itr++) {
      if (*big.itr != Xclude) {
        if (new_data.len + 1 > max_len)
          return false;
        new_data.push_back(*big.itr);
      }
    }

    if (*small.itr != Xclude && *small.itr != *big.itr) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(*small.itr);
    }
  }

  for (; big.itr != big.end; big.itr++) {
    if (*big.itr != Xclude) {
      if (new_data.len + 1 > max_len)
        return false;
      new_data.push_back(*big.itr);
    }
  }

  return true;
}

bool heuristic_merge(JData<vid_t> &new_data, size_t max_len,
    std::vector<SortedRange> &kid_itrs, vid_t Xclude)
{
  if (kid_itrs.size() < 32)
    return balance_line_merge(new_data, max_len, kid_itrs, Xclude);
  else
    return heap_merge(new_data, max_len, kid_itrs, Xclude);
}
