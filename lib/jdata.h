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

#include <algorithm>
#include <cstring>
#include "defs.h"

#ifndef NDEBUG
#define DEBUG_ONLY(statement) statement
#else
#define DEBUG_ONLY(statement)
#endif

template <typename DataType>
class JDataTable;

template <typename DataType, typename SizeType = esize_t>
struct JData {
  friend class JDataTable<DataType>;

private:
  inline JData(SizeType ml) : len(0) {
    DEBUG_ONLY(max_len = ml);
  }

public:
  DEBUG_ONLY(SizeType max_len);
  SizeType len; 
  DataType ele[1];
  #define SIZEOF_JDATA(DataType, len) (sizeof(JData<DataType>)-sizeof(DataType)+sizeof(DataType)*len)

  JData() = delete;
  ~JData() = default;

  JData(JData const &other) = delete;
  JData(JData &&other) = delete;
  JData& operator=(JData const &other) = delete;
  JData& operator=(JData &&other) = delete;

  inline void push_back(DataType X) {
    assert(len < max_len);
    ele[len++] = X;
  }

  inline SizeType size() const { return len; }

  inline DataType* begin() const { return const_cast<DataType*>(ele); }
  inline DataType* end() const { return const_cast<DataType*>(ele + len); }

  inline DataType const * cbegin() const { return ele; }
  inline DataType const * cend() const { return ele + len; }
  
  inline bool binary_search(DataType const X) const {
    return std::binary_search(begin(), end(), X);
  }

  inline bool binary_search(JData const &other) const {
    DataType *beg = begin();
    for (DataType const key : other) {
      beg = std::lower_bound(beg, end(), key);
      if (beg == end() || *beg != key)
        return false;
    }
    return true;
  }
};

template <typename DataType>
class JDataTable {
private:
  std::vector<size_t> offsets;

  size_t alloc_end;
  size_t alloc_max;
  char *allocation;

public:

  JDataTable() = delete;

  JDataTable(size_t max_offsets) :
    offsets(), alloc_end(0),
    alloc_max(SIZEOF_JDATA(DataType, 0) * max_offsets + sizeof(DataType) * max_offsets),
    allocation((char*)malloc(alloc_max))
  {
    if (alloc_max != 0) {
      if (allocation == nullptr) throw std::bad_alloc();
      offsets.reserve(max_offsets);
      new(allocation) JData<DataType>(0);
      alloc_end = SIZEOF_JDATA(DataType, 0);
    }
  }

  JDataTable(size_t max_offsets, size_t memory_limit) :
    offsets(), alloc_end(0),
    alloc_max(memory_limit), allocation((char*)malloc(alloc_max))
  {
    if (alloc_max != 0) {
      if (allocation == nullptr) throw std::bad_alloc();
      offsets.reserve(max_offsets);
      new(allocation) JData<DataType>(0);
      alloc_end = SIZEOF_JDATA(DataType, 0);
    }
  }

  JDataTable(JDataTable &&other) :
    offsets(std::move(other.offsets)), alloc_end(other.alloc_end),
    alloc_max(other.alloc_max), allocation(other.allocation)
  {
    other.alloc_end = 0;
    other.alloc_max = 0;
    other.allocation = nullptr;
  }

  JDataTable(JDataTable const &other) :
    offsets(other.offsets), alloc_end(other.alloc_end),
    alloc_max(other.alloc_max), allocation((char*)malloc(alloc_max))
  {
    if (alloc_max != 0) {
      if (allocation == nullptr) throw std::bad_alloc();
      offsets.reserve(other.offsets.capacity());
      memcpy(allocation, other.allocation, alloc_end);
    }
  }
  
  JDataTable(JDataTable const &other, size_t partial_end) :
    offsets(other.offsets.cbegin(), other.offsets.cbegin() + std::min(partial_end, other.size())),
    alloc_end(partial_end < other.size() ? other.offsets[partial_end] : 0),
    alloc_max(other.alloc_max), allocation((char*)malloc(alloc_max))
  {
    if (alloc_max != 0) {
      if (allocation == nullptr) throw std::bad_alloc();
      offsets.reserve(other.offsets.capacity());
      memcpy(allocation, other.allocation, alloc_end);
    }
  }

  JDataTable<DataType> & operator=(JDataTable<DataType> &&other)
  {
    offsets = std::move(other.offsets);

    alloc_end = other.alloc_end;
    other.alloc_end = 0;

    alloc_max = other.alloc_max;
    other.alloc_max = 0;

    free(allocation);
    allocation = other.allocation;
    other.allocation = nullptr;

    return *this;
  }

  JDataTable& operator=(JDataTable const &other) = delete;

  ~JDataTable() {
    free(allocation);
  }



  inline size_t size() const {
    return offsets.size();
  }

  inline size_t newJData(size_t max_len, bool const requires_max = true) {
    // If zero-length, use the dummy allocation (sparsifying).
    if (max_len == 0) {
      offsets.push_back(0);
      return size() - 1;
    }

    // Check remaining space.
    if (alloc_max - alloc_end < SIZEOF_JDATA(DataType, 0))
      throw std::bad_alloc();
    size_t const remaining_len =
      (alloc_max - alloc_end - SIZEOF_JDATA(DataType, 0)) / sizeof(DataType);
    if (max_len > remaining_len) {
      if (requires_max)
        throw std::bad_alloc();
      else
        max_len = remaining_len;
    }

    // Do the allocation.
    offsets.push_back(alloc_end);
    new(allocation + alloc_end) JData<DataType>(max_len);
    alloc_end += SIZEOF_JDATA(DataType, max_len);
    return size() - 1;
  }

  inline JData<DataType> & operator[](size_t const index) {
    assert(index < size());
    return *((JData<DataType> *)(allocation + offsets[index]));
  }

  inline JData<DataType> const & operator[](size_t const index) const {
    assert(index < size());
    return *((JData<DataType> *)(allocation + offsets[index]));
  }

  inline void shrinkJData(size_t const index) {
    assert(index == size() - 1);
    if (offsets[index] == 0) return;

    JData<DataType> &data = operator[](index);
    assert(offsets[index] + SIZEOF_JDATA(DataType, data.max_len) == alloc_end);
    if (data.len != 0) {
      DEBUG_ONLY(data.max_len = data.len);
      alloc_end = offsets[index] + SIZEOF_JDATA(DataType, data.len);
    } else {
      alloc_end = offsets[index];
      offsets[index] = 0;
    }
  }
  
  inline void deleteJData(size_t const index) {
    assert(index == size() - 1);
    if (offsets[index] != 0) {
      JData<DataType> &data = operator[](index);
      assert(offsets[index] + SIZEOF_JDATA(DataType, data.max_len) == alloc_end);
      alloc_end = offsets[index];
    }
    offsets.pop_back();
  }
};

#undef SIZEOF_JDATA
#undef DEBUG_ONLY
