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

bool balance_line_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
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
bool heap_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
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

bool asymmetric_merge(JData<vid_t> &new_data, std::vector<SortedRange> &kid_itrs,
    vid_t Xclude, size_t max_len)
{
  assert(kid_itrs.size() == 2);
  SortedRange big = kid_itrs[0];
  SortedRange small = kid_itrs[1];
  if (big.size() < small.size())
    std::swap(big,small);
  if (big.size() < small.size() * 8)
    return balance_line_merge(new_data, kid_itrs, Xclude, max_len);

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

