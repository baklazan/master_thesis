#ifndef MATRIX_BAND_H
#define MATRIX_BAND_H
#include <vector>
#include <algorithm>
#include <cstdio> //stderr

template<class T>
class MatrixBand {
private:
  std::vector<int> row_starts_;
  std::vector<int> row_ends_;
  T out_of_bounds_value_backup_;
  T out_of_bounds_value_;
  int matrix_dimension_;
  std::vector<std::vector<T>> values_;
public:
  MatrixBand(std::vector<int> row_starts, std::vector<int> row_ends, T default_value, T out_of_bounds_value_)
    : row_starts_(row_starts), row_ends_(row_ends), out_of_bounds_value_backup_(out_of_bounds_value_), 
      matrix_dimension_(static_cast<int>(row_starts.size())), values_(matrix_dimension_) {
    for (int y = 0; y < matrix_dimension_; y++) {
      int row_length = row_ends_[y] - row_starts_[y];
      values_[y].resize(row_length, default_value);
    }
  }
  
  T &Get(int row, int column) {
    if (row < 0 || row >= matrix_dimension_ || column < row_starts_[row] || column >= row_ends_[row]) {
      out_of_bounds_value_ = out_of_bounds_value_backup_;
      fprintf(stderr, "out of bounds access at %d %d\n", row, column);
      return out_of_bounds_value_;
    }
    else {
      return values_[row][column - row_starts_[row]];
    }
  }
  
  class MatrixRow {
    friend class MatrixBand;
  private:
    int row_;
    MatrixBand *parent_matrix_;
    MatrixRow(int row, MatrixBand *parent_matrix) : row_(row), parent_matrix_(parent_matrix) {}
  public:
    T &operator[](int index) {
      return parent_matrix_->Get(row_, index);
    }
  };
  
  MatrixRow operator[](int index) {
    return MatrixRow(index, this);
  }
  
  T &operator[](std::pair<int, int> indices) {
    return Get(indices.first, indices.second);
  }
};

#endif
