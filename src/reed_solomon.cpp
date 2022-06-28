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

#include "reed_solomon.h"

#include <array>
#include <bitset>
#include <cassert>
#include <limits>

#include "cpuid.hpp"

// Initialize static member data
const CpuId::CpuIdInternal CpuId::cpu_info_;
ReedSolomon::Galois ReedSolomon::galois_;

namespace {
size_t GetShardSize(const std::vector<std::vector<uint8_t>>& shards) {
  for (const auto& e : shards) {
    if (0 != e.size()) {
      return e.size();
    }
  }
  return 0;
}

size_t GetShardSize(
    const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards) {
  for (const auto& e : shards) {
    if (0 != e.get().size()) {
      return e.get().size();
    }
  }
  return 0;
}

bool CheckShards(const std::vector<std::vector<uint8_t>>& shards,
                 std::size_t rows,
                 bool allow_empty) {
  std::size_t shard_size = GetShardSize(shards);
  assert(0 != shard_size);
  if (0 == shard_size) {
    return false;
  }
  std::size_t r = 0;
  for (const auto& e : shards) {
    ++r;
    if (rows <= r) {
      break;
    }
    if (e.size() != shard_size) {
      if (0 != e.size() || !allow_empty) {
        return false;
      }
    }
  }
  return true;
}

bool CheckShards(
    const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards,
    std::size_t rows,
    bool allow_empty) {
  std::size_t shard_size = GetShardSize(shards);
  assert(0 != shard_size);
  if (0 == shard_size) {
    return false;
  }
  std::size_t r = 0;
  for (const auto& e : shards) {
    ++r;
    if (rows <= r) {
      break;
    }
    if (e.get().size() != shard_size) {
      if (0 != e.get().size() || !allow_empty) {
        return false;
      }
    }
  }
  return true;
}
}  // namespace

ReedSolomon::ReedSolomon(std::size_t k, std::size_t r, MatrixType matrix_type)
    : k_(k), r_(r) {
  // Once initialized, exp_table_[0] is 1
  assert(1 == galois_.exp_table_[0]);

  assert(0 < k && k <= std::numeric_limits<gf>::max());
  assert(0 < r && r <= std::numeric_limits<gf>::max());

  if (k < 1 || (std::numeric_limits<gf>::max)() < k) {
    throw std::out_of_range("Invalid data shards.");
  }
  if (r < 1 || std::numeric_limits<gf>::max() < r) {
    throw std::out_of_range("Invalid parity shards.");
  }
  std::size_t shards = k + r;
  assert(shards <= std::numeric_limits<gf>::max());
  if (std::numeric_limits<gf>::max() < shards) {
    throw std::out_of_range("Invalid total shards.");
  }

  switch (matrix_type) {
    case ReedSolomon::MatrixType::kCauchy:
      matrix_ = BuildMatrixCauchy(shards, k);
      break;
    case ReedSolomon::MatrixType::kParV1:
      if (1 == r) {
        matrix_ = BuildMatrixVandermondeFast(shards, k);
      } else {
        matrix_ = BuildMatrixParV1(shards, k);
      }
      break;
    case ReedSolomon::MatrixType::kVandermonde:
      [[fallthrough]];
    default:
      if (1 == r) {
        matrix_ = BuildMatrixVandermondeFast(shards, k);
      } else {
        matrix_ = BuildMatrixVandermonde(shards, k);
      }
      break;
  }

#if _DEBUG
  std::cout << "Matrix:\n";
  for (std::size_t r = 0; r < matrix_.Rows(); ++r) {
    for (std::size_t c = 0; c < matrix_.Columns(); ++c) {
      std::cout << static_cast<int>(matrix_.At(r, c)) << " ";
    }
    std::cout << '\n';
  }
#endif

  if (CpuId::AVX2()) {
  }
}

bool ReedSolomon::Encode(std::vector<std::vector<uint8_t>>& shards) {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }
  if (!CheckShards(shards, GetDataShardCount(), false)) {
    assert(false);
    return false;
  }

  std::vector<std::reference_wrapper<std::vector<uint8_t>>> input(
      shards.begin(), shards.begin() + GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> output(
      shards.begin() + GetDataShardCount(), shards.end());
  CodeSomeShards(matrix_, input, output);
  return true;
}

bool ReedSolomon::Encode(
    std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards) {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }
  if (!CheckShards(shards, GetDataShardCount(), false)) {
    assert(false);
    return false;
  }

  std::vector<std::reference_wrapper<std::vector<uint8_t>>> input(
      shards.begin(), shards.begin() + GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> output(
      shards.begin() + GetDataShardCount(), shards.end());
  CodeSomeShards(matrix_, input, output);
  return true;
}

bool ReedSolomon::Reconstruct(std::vector<std::vector<uint8_t>>& shards) {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }

  std::size_t shard_size = GetShardSize(shards);
  assert(0 != shard_size);
  if (0 == shard_size) {
    return false;
  }
  std::size_t missing = 0;
  for (std::size_t d = 0; d < GetDataShardCount(); ++d) {
    if (shards[d].size() != shard_size) {
      assert(0 == shards[d].size());
      if (0 != shards[d].size()) {
        return false;
      }
      ++missing;
    }
  }
  if (0 == missing) {
    return true;
  }
  if (GetParityShardCount() < missing) {
    assert(false);
    return false;
  }

  Matrix present_matrix(GetDataShardCount(), GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> present_shards;

  present_shards.reserve(GetDataShardCount());
  std::size_t present = 0;
  for (std::size_t r = 0;
       r < GetTotalShardCount() && present < GetDataShardCount(); ++r) {
    if (shards[r].size() == shard_size) {
      present_matrix.SetRow(present, matrix_.GetRow(r));
      present_shards.emplace_back(std::ref(shards[r]));
      ++present;
    }
  }
  Matrix decode_matrix = present_matrix.Invert();
  // TO-DO: cache decode_matrix, key is the missing indices

  // Re-create
  Matrix missing_matrix(GetTotalShardCount(), GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> missing_rows;
  missing_rows.reserve(missing);
  missing = 0;
  for (std::size_t r = 0; r < GetDataShardCount(); ++r) {
    if (shards[r].empty()) {
      missing_matrix.SetRow(present + missing, decode_matrix.GetRow(r));
      missing_rows.emplace_back(std::ref(shards[r]));
      ++missing;
    }
  }
  CodeSomeShards(missing_matrix, present_shards, missing_rows);

  return true;
}

bool ReedSolomon::Reconstruct(
    std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards) {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }

  std::size_t shard_size = GetShardSize(shards);
  assert(0 != shard_size);
  if (0 == shard_size) {
    return false;
  }
  std::size_t missing = 0;
  for (std::size_t d = 0; d < GetDataShardCount(); ++d) {
    if (shards[d].get().size() != shard_size) {
      assert(0 == shards[d].get().size());
      if (0 != shards[d].get().size()) {
        return false;
      }
      ++missing;
    }
  }
  if (0 == missing) {
    return true;
  }
  if (GetParityShardCount() < missing) {
    assert(false);
    return false;
  }

  Matrix present_matrix(GetDataShardCount(), GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> present_shards;
  present_shards.reserve(GetDataShardCount());
  std::size_t present = 0;
  for (std::size_t r = 0;
       r < GetTotalShardCount() && present < GetDataShardCount(); ++r) {
    if (shards[r].get().size() == shard_size) {
      present_matrix.SetRow(present, matrix_.GetRow(r));
      present_shards.emplace_back(std::ref(shards[r]));
      ++present;
    }
  }
  Matrix decode_matrix = present_matrix.Invert();
  // TO-DO: cache decode_matrix, key is the missing indices

  // Re-create
  Matrix missing_matrix(GetTotalShardCount(), GetDataShardCount());
  std::vector<std::reference_wrapper<std::vector<uint8_t>>> missing_rows;
  missing_rows.reserve(missing);
  missing = 0;
  for (std::size_t r = 0; r < GetDataShardCount(); ++r) {
    if (0 == shards[r].get().size()) {
      missing_matrix.SetRow(present + missing, decode_matrix.GetRow(r));
      missing_rows.emplace_back(std::ref(shards[r]));
      ++missing;
    }
  }
  CodeSomeShards(missing_matrix, present_shards, missing_rows);

  return true;
}

bool ReedSolomon::Verify(
    const std::vector<std::vector<uint8_t>>& shards) const noexcept {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }
  if (!CheckShards(shards, GetTotalShardCount(), false)) {
    return false;
  }
  // TO-DO
  assert(false);
  return true;
}

bool ReedSolomon::Verify(
    const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& shards)
    const noexcept {
  assert(shards.size() == GetTotalShardCount());
  if (shards.size() != GetTotalShardCount()) {
    return false;
  }
  if (!CheckShards(shards, GetTotalShardCount(), false)) {
    return false;
  }
  // TO-DO
  assert(false);
  return true;
}

// private
ReedSolomon::Matrix ReedSolomon::Vandermonde(std::size_t rows,
                                             std::size_t columns) noexcept {
  assert(columns < rows);
  assert(rows <= std::numeric_limits<gf>::max());

  Matrix matrix(rows, columns);
  for (gf r = 0; r < rows; ++r) {
    for (gf c = 0; c < columns; ++c) {
      matrix.At(r, c) = galois_.Exp(r, c);
    }
  }
  return matrix;
}

ReedSolomon::Matrix ReedSolomon::BuildMatrixCauchy(
    std::size_t rows,
    std::size_t columns) noexcept {
  assert(columns < rows);
  assert(rows <= std::numeric_limits<gf>::max());

  Matrix matrix(rows, columns);
  for (gf r = 0; r < columns; ++r) {
    matrix.At(r, r) = 1;
  }
  for (gf r = static_cast<gf>(columns); r < rows; ++r) {
    for (gf c = 0; c < columns; ++c) {
      matrix.At(r, c) = galois_.inverse_table_[galois_.Subtract(r, c)];
    }
  }
  return matrix;
}

ReedSolomon::Matrix ReedSolomon::BuildMatrixParV1(
    std::size_t rows,
    std::size_t columns) noexcept {
  assert(columns < rows);
  assert(rows <= std::numeric_limits<gf>::max());

  Matrix matrix(rows, columns);
  for (gf r = 0; r < columns; ++r) {
    matrix.At(r, r) = 1;
  }
  for (gf r = static_cast<gf>(columns); r < rows; ++r) {
    for (gf c = 0; c < columns; ++c) {
      matrix.At(r, c) = galois_.Exp(c + 1, r - static_cast<gf>(columns));
    }
  }
  return matrix;
}

// When r == 1
ReedSolomon::Matrix ReedSolomon::BuildMatrixVandermondeFast(
    std::size_t rows,
    std::size_t columns) noexcept {
  assert(columns + 1 == rows);

  Matrix matrix(rows, columns);
  for (std::size_t r = 0; r < columns; ++r) {
    matrix.At(r, r) = 1;
  }
  for (std::size_t c = 0; c < columns; ++c) {
    matrix.At(columns, c) = 1;
  }
  return matrix;
}

ReedSolomon::Matrix ReedSolomon::BuildMatrixVandermonde(
    std::size_t rows,
    std::size_t columns) noexcept(false) {
  assert(columns < rows);

  auto vm = Vandermonde(rows, columns);
  auto top = vm.Sub(0, 0, columns, columns);
  return vm.Multiply(top.Invert());
}

void ReedSolomon::CodeSomeShards(
    Matrix matrix,
    const std::vector<std::reference_wrapper<std::vector<uint8_t>>>& input,
    const std::vector<std::reference_wrapper<std::vector<uint8_t>>>&
        output) noexcept {
  assert(input.size() > 0);
  assert(output.size() > 0);
  // First row -> output rows
  // Note: [] is fast then at()
  {
    for (std::size_t output_index = 0; output_index < output.size();
         ++output_index) {
      auto& output_raw = output[output_index].get();
      auto& input_raw = input[0].get();
      if (output_raw.size() != input_raw.size()) {
        output_raw.resize(input_raw.size());
      }
      for (std::size_t c = 0; c < input_raw.size(); ++c) {
        output_raw[c] = galois_.Multiply(
            matrix.At(input.size() + output_index, 0), input_raw[c]);
      }
    }
  }

  for (std::size_t input_index = 1; input_index < input.size(); ++input_index) {
    auto& input_raw = input[input_index].get();
    for (std::size_t output_index = 0; output_index < output.size();
         ++output_index) {
      auto& output_raw = output[output_index].get();
      for (std::size_t c = 0; c < output_raw.size(); ++c) {
        output_raw[c] =
            galois_.Add(output_raw[c],
                        galois_.Multiply(
                            matrix.At(input.size() + output_index, input_index),
                            input_raw[c]));
      }
    }
  }
}

/*
 * The polynomial used to generate the logarithm table.
 *
 * Markdown: $$...$$ or ```math...``` or `$...$`
 *
 * $$
 * f(x) = \alpha_nx^n + \alpha_{n-1}x^{n-1} + ... + \alpha_1x + \alpha_0
 * \alpha_i \in 0, 1 (i = 0, ..., n)
 * $$
 *
 * \alpha elements map to 0,...,2^n-1
 *
 * Primitive polynomial (x -> D):
 * \alpha is the primitive element of GF(2^n)
 *
 * $$
 * f(D) = D^n + \alpha_{n-1}D^{n-1} + ... + \alpha_1D + 1
 * $$
 *
 * https://octave-online.net/
 *
 * octave:1> pkg load communications; primpoly(8,'all');
 * octave:2> for i = 2:16 primpoly(i); end
 *
 * For 8-bit Galois Field, we use 0b0001'1101 (29, 0x1d)
 * The possibilities see octave primpoly(8,'all')
 */
constexpr std::array<std::bitset<8>, 17> kPrimitivePolynomial = {
    0, 0,
    // bitset<8> is just enough. The highest and the lowest bit(\alpha_0) are
    // always 1, so the highest bit is ignored here.
    0b0000'0011,  // 2  D^2+D+1
    0b0000'0011,  // 3  D^3+D+1
    0b0000'0011,  // 4  D^4+D+1
    0b0000'0101,  // 5  D^5+D^2+1
    0b0000'0011,  // 6  D^6+D+1
    0b0000'0011,  // 7  D^7+D+1
    0b0001'1101,  // 8  D^8+D^4+D^3+D^2+1
    0b0001'0001,  // 9  D^9+D^4+1
    0b0000'1001,  // 10 D^10+D^3+1
    0b0000'0101,  // 11 D^11+D^2+1
    0b0101'0011,  // 12 D^12+D^6+D^4+D+1
    0b0001'1011,  // 13 D^13+D^4+D^3+D+1
    0b0010'1011,  // 14 D^14+D^5+D^3+D+1
    0b0000'0011,  // 15 D^15+D+1
    0b0010'1101,  // 16 D^16+D^5+D^3+D^2+1

    /*
     * Other implementation use:
     * 7:  D^7+D^3+1
     * 14: D^14+D^10+D^6+D+1
     * 16: D^16+D^12+D^3+D+1
     */
};

ReedSolomon::Galois::Galois() {
  // ReedSolomon::galois_ is static, so this runs once.
  GenerateGfTables();
  GenerateMultiplicationTable();
}

void ReedSolomon::Galois::GenerateGfTables() noexcept {
  gf mask = 1;  // \alpha_0D^0 = 1
  exp_table_[kGfBits] = 0;
  for (std::size_t i = 0; i < kGfBits; ++i, mask <<= 1) {
    exp_table_[i] = mask;  // exp_table_[i] = D^i
    // log is the inverse operation of exp
    log_table_[exp_table_[i]] = static_cast<uint32_t>(i);

    if (kPrimitivePolynomial.at(kGfBits)[i]) {
      exp_table_[kGfBits] ^= mask;
    }
  }

  // exp_table_[kGfBits] = D^kGfBits is complete, compute its inverse.
  log_table_[exp_table_[kGfBits]] = kGfBits;

  mask = 1 << (kGfBits - 1);
  for (std::size_t i = kGfBits + 1; i < kGfCycle; ++i) {
    if (exp_table_[i - 1] >= mask) {
      exp_table_[i] = exp_table_[kGfBits] ^ ((exp_table_[i - 1] ^ mask) << 1);
    } else {
      exp_table_[i] = exp_table_[i - 1] << 1;
    }
    log_table_[exp_table_[i]] = static_cast<uint32_t>(i);
  }

  // For efficiency, exp_table_[] has 2 copies, so that a simple multiplication
  // of two numbers can be resolved without mod(%) operation.
  for (std::size_t i = 0; i < kGfCycle; ++i) {
    exp_table_[i + kGfCycle] = exp_table_[i];
  }

  // log_table_[0] and inverse_table_[0] is undefined
  for (std::size_t i = 1; i < kFieldSize; ++i) {
    // inverse_table_[\alpha^n] = \alpha^{kGfCycle - n}
    // i = \alpha^n, so n = log(i)
    inverse_table_[i] = exp_table_[kGfCycle - log_table_[i]];
  }

  // Unit Tests
  if constexpr (8 == kGfBits) {
    // The same as https://github.com/templexxx/reedsolomon
    assert(0x01 == exp_table_[0x00] && 0x01 == exp_table_[0xff]);
    assert(0x8e == exp_table_[0xfe] && 0x8e == exp_table_[0x1fd]);

    assert(0x00 == log_table_[0x01]);
    assert(0xaf == log_table_[0xff]);

    assert(0x01 == inverse_table_[0x01]);
    assert(0xfd == inverse_table_[0xff]);
  }
}

void ReedSolomon::Galois::GenerateMultiplicationTable() noexcept {
  auto mod = [](int x) -> gf {
    while (x >= kGfCycle) {
      x -= kGfCycle;
      x = (x >> kGfBits) + (x & kGfCycle);
    }
    return static_cast<gf>(x);
  };

  // log_table_[0] is undefined
  for (std::size_t i = 1; i < kFieldSize; ++i) {
    for (std::size_t j = 1; j < kFieldSize; ++j) {
      // warning C28020 is known
      mul_table_[(i << kGfBits) | j] =
          exp_table_[mod(log_table_[i] + log_table_[j])];
    }
  }

  // Unit Tests
  if constexpr (8 == kGfBits) {
    // The same as https://github.com/templexxx/reedsolomon
    assert(0x01 == mul_table_[0x0101]);
    assert(0xff == mul_table_[0x01ff]);
    assert(0xff == mul_table_[0xff01]);
    assert(0xe2 == mul_table_[0xffff]);
  }
}

ReedSolomon::gf ReedSolomon::Galois::Exp(gf x, gf n) const noexcept {
  // $$0^0 == 1$$
  if (0 == n) {
    return 1;
  }
  if (0 == x) {
    return 0;
  }

  std::size_t result = log_table_[x] * static_cast<size_t>(n);
  while (exp_table_.size() <= result) {
    result -= exp_table_.size();
  }
  // warning C28020 is known
  return exp_table_[result];
}

ReedSolomon::Matrix ReedSolomon::Matrix::Augment(
    const Matrix& right) const noexcept {
  assert(Rows() == right.Rows());

  Matrix result(Rows(), Columns() + right.Columns());
  for (std::size_t r = 0; r < Rows(); ++r) {
    for (std::size_t c = 0; c < Columns(); c++) {
      result.At(r, c) = At(r, c);
    }
    for (std::size_t c = 0; c < right.Columns(); c++) {
      result.At(r, Columns() + c) = right.At(r, c);
    }
  }
  return result;
}

ReedSolomon::Matrix ReedSolomon::Matrix::Sub(std::size_t rmin,
                                             std::size_t cmin,
                                             std::size_t rmax,
                                             std::size_t cmax) const noexcept {
  assert(rmin < rmax && rmax - rmin <= Rows());
  assert(cmin < cmax && cmax - cmin <= Columns());

  Matrix result(rmax - rmin, cmax - cmin);
  for (std::size_t r = rmin; r < rmax; ++r) {
    for (std::size_t c = cmin; c < cmax; ++c) {
      result.At(r - rmin, c - cmin) = At(r, c);
    }
  }
  return result;
}

ReedSolomon::Matrix ReedSolomon::Matrix::Multiply(
    const Matrix& right) const noexcept {
  assert(Columns() == right.Rows());

  Matrix result(Rows(), right.Columns());
  for (std::size_t r = 0; r < Rows(); ++r) {
    for (std::size_t c = 0; c < right.Columns(); ++c) {
      gf value = 0;
      for (std::size_t i = 0; i < Columns(); i++) {
        value = galois_.Add(value, galois_.Multiply(At(r, i), right.At(i, c)));
      }
      result.At(r, c) = value;
    }
  }
  return result;
}

ReedSolomon::Matrix ReedSolomon::Matrix::Invert() const noexcept(false) {
  // Only square matrices can be inverted
  assert(IsSquare());
  if (!IsSquare()) {
    throw std::logic_error("Matrix is not square");
  }

  // Create a matrix by augmenting this one with an identity matrix on the
  // right.
  Matrix matrix = Augment(Identity(Rows()));
  // Do Gauss-Jordan elimination to transform the left half into
  // an identity matrix.
  if (!matrix.GaussianElimination()) {
    throw std::logic_error("Matrix is singular");
  }
  // The right half is now the inverse.
  return matrix.Sub(0, Columns(), Rows(), Columns() * 2);
}

/* Another solution: Gauss-Jordan Elimination.
 * https://www.avrfreaks.net/sites/default/files/forum_attachments/Gauss-Jordan%20code.pdf
 * However, its principal weaknesses are (i) that it requires all the right-hand
 * sides to be stored and manipulated at the same time, and (ii) that when the
 * inverse matrix is not desired, Gauss-Jordan is three times slower than the
 * best alternative technique for solving a single linear set (¡ì2.3). The
 * method¡¯s principal strength is that it is as stable as any other direct
 * method, perhaps even a bit more stable when full pivoting is used (see
 * below).
 *
 * 1. upper triangular matrix. $$A\vec{x} = \vec{b} \to U\vec{x} = \vec{c}$$
 * 2. back substitution
 *
 * Return false if singular, which means det(A) = 0.
 */
bool ReedSolomon::Matrix::GaussianElimination() noexcept {
  // Clear out the part below the main diagonal and scale the main diagonal to
  // be 1.
  for (std::size_t r = 0; r < Rows(); ++r) {
    // If the element on the diagonal is 0, find a row below that has a non-zero
    // and swap them.
    if (0 == At(r, r)) {
      for (std::size_t row_below = r + 1; row_below < Rows(); ++row_below) {
        if (0 != At(row_below, r)) {
          SwapRows(r, row_below);
          break;
        }
      }
    }
    // If we couldn't find one, the matrix is singular.
    if (0 == At(r, r)) {
      return false;
    }
    // Scale to 1.
    if (1 != At(r, r)) {
      gf scale = galois_.Divide(1, At(r, r));
      for (std::size_t c = 0; c < Columns(); ++c) {
        At(r, c) = galois_.Multiply(At(r, c), scale);
      }
    }
    // Make everything below the 1 be a 0 by subtracting a multiple of it.
    for (std::size_t row_below = r + 1; row_below < Rows(); ++row_below) {
      gf scale = At(row_below, r);
      if (0 != scale) {
        for (std::size_t c = 0; c < Columns(); ++c) {
          At(row_below, c) = galois_.Subtract(
              At(row_below, c), galois_.Multiply(At(r, c), scale));
        }
      }
    }
  }

  // Now clear the part above the main diagonal.
  for (std::size_t d = 1; d < Rows(); ++d) {
    for (std::size_t row_above = 0; row_above < d; ++row_above) {
      gf scale = At(row_above, d);
      if (0 != scale) {
        for (std::size_t c = 0; c < Columns(); ++c) {
          At(row_above, c) = galois_.Subtract(
              At(row_above, c), galois_.Multiply(At(d, c), scale));
        }
      }
    }
  }
  return true;
}
