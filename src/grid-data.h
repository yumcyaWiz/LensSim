#ifndef _GRID_DATA_H
#define _GRID_DATA_H
#include "core/type.h"

using namespace Prl2;

template <typename T>
class GridData {
 public:
  const unsigned int nrows;
  const unsigned int ncols;

  T* data;

  GridData(unsigned int _nrows, unsigned int _ncols)
      : nrows(_nrows), ncols(_ncols) {
    data = new T[ncols * nrows];
  }
  ~GridData() { delete[] data; }

  T get(unsigned int i, unsigned int j) const { return data[j + ncols * i]; }

  void set(unsigned int i, unsigned int j, const T& value) {
    data[j + ncols * i] = value;
  }
};
#endif