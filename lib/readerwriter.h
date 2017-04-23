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

#include <fstream>

#include "defs.h"

struct xs1 {
  unsigned tail;
  unsigned head;
  float weight;
};

class XS1Reader {
private:
  std::ifstream stream;
  xs1 buf;
public:
  XS1Reader(char const *const filename) :
    stream(filename, std::ios::binary) {}

  bool read(vid_t &X, vid_t &Y) {
    if (!stream.eof()) {
      stream.read((char*)&buf, sizeof(xs1));
      X = buf.tail;
      Y = buf.head;
      return true;
    }
    return false;
  }
};

class XS1Writer {
private:
  std::ofstream stream;
  xs1 buf;
public:
  XS1Writer(char const *const filename) :
    stream(filename, std::ios::binary | std::ios::trunc) {
    buf.weight = 1.0;
  }

  void write(vid_t const X, vid_t const Y) {
    buf.tail = X;
    buf.head = Y;
    stream.write((char*)&buf, sizeof(xs1));
  }
};

class SNAPReader {
private:
  std::ifstream stream;
public:
  SNAPReader(char const *const filename) :
    stream(filename) {}

  bool read(vid_t &X, vid_t &Y) {
    bool result = (bool) (stream >> X);
    result &= (bool) (stream >> Y);
    return result;
  }
};

class SNAPWriter {
private:
  std::ofstream stream;
public:
  SNAPWriter(char const *const filename) :
    stream(filename, std::ios::trunc) {}

  void write(vid_t const X, vid_t const Y) {
    stream << X << ' ' << Y << std::endl;
  }
};
