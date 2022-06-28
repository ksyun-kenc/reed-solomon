/*
 * Copyright 2021-present Ksyun
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

class ReedSolomon {
 public:
  /*
   * The following parameter defines how many bits are used for field elements.
   * The code supports any value from 2 to 16 but fastest operation is achieved
   * with 8 bit elements.
   * This is the only parameter you may want to change.
   */
  static const std::size_t kGfBits = 8;

  // You should not need to change anything beyond this point.
  using gf = std::conditional_t<kGfBits <= 8, std::uint8_t, std::uint16_t>;

  // default: kVandermonde
  enum class MatrixType { kVandermonde, kCauchy, kParV1 };

  // k: count of data shards
  // r: count of parity shards
  ReedSolomon(std::size_t k,
              std::size_t r,
              MatrixType matrix_type = MatrixType::kVandermonde);
  ReedSolomon(const ReedSolomon& right) noexcept = default;
  ReedSolomon(ReedSolomon&& right) noexcept = default;
  ~ReedSolomon() {}
  ReedSolomon& operator=(const ReedSolomon& right) noexcept {
    if (this != std::addressof(right)) {
      k_ = right.k_;
      r_ = right.r_;
      matrix_ = right.matrix_;
    }
    return *this;
  }
  ReedSolomon& operator=(ReedSolomon&& right) noexcept {
    if (this != std::addressof(right)) {
      k_ = right.k_;
      r_ = right.r_;
      matrix_ = std::move(right.matrix_);
    }
    return *this;
  }

  bool Encode(std::vector<std::vector<uint8_t>>& shards);
  bool Encode(
      std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards);

  bool Reconstruct(std::vector<std::vector<uint8_t>>& shards);
  bool Reconstruct(
      std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards);

  bool Verify(const std::vector<std::vector<uint8_t>>& shards) const noexcept;
  bool Verify(const std::vector<std::reference_wrapper<std::vector<uint8_t>>>&
                  shards) const noexcept;

  std::size_t GetDataShardCount() const noexcept { return k_; }
  std::size_t GetParityShardCount() const noexcept { return r_; }
  std::size_t GetTotalShardCount() const noexcept { return matrix_.Rows(); }

 private:
  class Matrix {
   public:
    Matrix() noexcept : rows_(0), columns_(0) {}

    Matrix(std::size_t rows, std::size_t columns) noexcept
        : rows_(rows), columns_(columns), data_(rows * columns) {}

    Matrix(std::size_t rows,
           std::size_t columns,
           std::initializer_list<gf> list) noexcept
        : rows_(rows), columns_(columns), data_{std::move(list)} {
      if (rows * columns != data_.size()) {
        data_.resize(rows * columns);
      }
    }

    Matrix(const Matrix& right) noexcept
        : rows_(right.rows_), columns_(right.columns_), data_(right.data_) {}

    Matrix(Matrix&& right) noexcept
        : rows_(right.rows_),
          columns_(right.columns_),
          data_(std::move(right.data_)) {}

    Matrix& operator=(const Matrix& right) noexcept {
      if (this != std::addressof(right)) {
        rows_ = right.rows_;
        columns_ = right.columns_;
        data_ = right.data_;
      }
      return *this;
    }
    Matrix& operator=(Matrix&& right) noexcept {
      if (this != std::addressof(right)) {
        rows_ = right.rows_;
        columns_ = right.columns_;
        data_ = std::move(right.data_);
      }
      return *this;
    }

    static Matrix Identity(std::size_t n) noexcept {
      Matrix matrix(n, n);
      for (std::size_t i = 0; i < n; ++i) {
        matrix.At(i, i) = 1;
      }
      return matrix;
    }

    const std::vector<gf>& Data() const noexcept { return data_; }
    std::vector<gf>& Data() noexcept { return data_; }
    std::size_t Rows() const noexcept { return rows_; }
    std::size_t Columns() const noexcept { return columns_; }
    bool IsSquare() const noexcept { return Rows() == Columns(); }

    gf& At(std::size_t r, std::size_t c) noexcept {
      assert(0 <= r && r < Rows());
      assert(0 <= c && c < Columns());
      return data_[r * Columns() + c];
    }
    const gf& At(std::size_t r, std::size_t c) const noexcept {
      assert(0 <= r && r < Rows());
      assert(0 <= c && c < Columns());
      return data_[r * Columns() + c];
    }

    std::vector<gf> GetRow(std::size_t r) noexcept {
      assert(0 <= r && r < Rows());
      std::vector<gf> row(data_.cbegin() + r * Columns(),
                          data_.cbegin() + (r + 1) * Columns());
      return row;
    }
    void SetRow(std::size_t r, const std::vector<gf> row) noexcept {
      assert(0 <= r && r < Rows());
      assert(Columns() == row.size());
      std::copy(row.cbegin(), row.cend(), data_.begin() + r * Columns());
    }

    void SwapRows(std::size_t r1, std::size_t r2) noexcept {
      assert(0 <= r1 && r1 < Rows());
      assert(0 <= r2 && r2 < Rows());
      for (std::size_t c = 0; c < Columns(); ++c) {
        std::swap(At(r1, c), At(r2, c));
      }
    }

    Matrix Augment(const Matrix& right) const noexcept;
    Matrix Sub(std::size_t rmin,
               std::size_t cmin,
               std::size_t rmax,
               std::size_t cmax) const noexcept;
    Matrix Multiply(const Matrix& right) const noexcept;
    Matrix Invert() const noexcept(false);
    bool GaussianElimination() noexcept;

   private:
    std::size_t rows_;
    std::size_t columns_;
    std::vector<gf> data_;
  };

  class Galois {
   public:
    Galois();
    void GenerateGfTables() noexcept;
    void GenerateMultiplicationTable() noexcept;
    // Subtraction and addition are both exclusive or in the Galois field.
    gf Add(gf a, gf b) const noexcept { return a ^ b; }
    gf Subtract(gf a, gf b) const noexcept { return a ^ b; }
    gf Multiply(gf a, gf b) const noexcept {
      if (0 == a || 0 == b) {
        return 0;
      }
      // warning C26451 is known
      return exp_table_[log_table_[a] + log_table_[b]];
    }
    gf Divide(gf a, gf b) const noexcept {
      if (0 == a) {
        return 0;
      }
      if (a == b) {
        return 1;
      }
      assert(0 != b);
      if (log_table_[a] < log_table_[b]) {
        // log_table_[0] is undefined, kGfCycle == log_table_.size() - 1
        // warning C26451 is known
        return exp_table_[kGfCycle - log_table_[b] + log_table_[a]];
      } else {
        return exp_table_[log_table_[a] - log_table_[b]];
      }
    }
    gf Exp(gf x, gf n) const noexcept;  // $$x^n$$

   public:
    static const std::size_t kFieldSize = 1 << kGfBits;  // count of elements
    static const gf kGfCycle = kFieldSize - 1;  // $$g^k=g^(k%kGfCycle)$$

    std::array<gf, 2 * kGfCycle> exp_table_;
    std::array<uint32_t, kFieldSize> log_table_;  // log_table_[0] is undefined
    std::array<gf, kFieldSize>
        inverse_table_;  // inverse_table_[0] is undefined
    std::array<gf, kFieldSize * kFieldSize> mul_table_;
  };

  static Matrix Vandermonde(std::size_t rows, std::size_t columns) noexcept;

  static Matrix BuildMatrixCauchy(std::size_t rows,
                                  std::size_t columns) noexcept;
  static Matrix BuildMatrixParV1(std::size_t rows,
                                 std::size_t columns) noexcept;
  static Matrix BuildMatrixVandermondeFast(std::size_t rows,
                                           std::size_t columns) noexcept;
  static Matrix BuildMatrixVandermonde(std::size_t rows,
                                       std::size_t columns) noexcept(false);

  static void CodeSomeShards(
      Matrix matrix,
      const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& input,
      const std::vector<std::reference_wrapper<std::vector<uint8_t>>>&
          output) noexcept;

 private:
  static Galois galois_;

  std::size_t k_;
  std::size_t r_;
  Matrix matrix_;
};
