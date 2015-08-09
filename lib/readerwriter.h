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
    bool result = (stream >> X);
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
