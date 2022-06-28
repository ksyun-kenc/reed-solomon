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

#include <algorithm>
#include <iostream>

#include "reed_solomon.h"

void Print(const std::vector<uint8_t>& shard) {
  std::cout << "  [";
  for (size_t c = 0; c < shard.size(); ++c) {
    std::cout << static_cast<int>(shard.at(c));
    if (c < shard.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";
}

void Print(const std::vector<std::vector<uint8_t>>& shards) {
  std::cout << "[\n";
  for (size_t r = 0; r < shards.size(); ++r) {
    Print(shards.at(r));
  }
  std::cout << "]\n";
}

int main() {
  const std::size_t kDataShardCount = 3;
  const std::size_t kParityShardCount = 1;
  const std::size_t kTotalShardCount = kDataShardCount + kParityShardCount;
  const std::size_t kShardSize = 64;

  // ReedSolomon rs(kDataShardCount, kParityShardCount,
  // ReedSolomon::MatrixType::kCauchy);
  ReedSolomon rs(kDataShardCount, kParityShardCount,
                 ReedSolomon::MatrixType::kVandermonde);
  std::vector<std::vector<uint8_t>> shards(kTotalShardCount,
                                           std::vector<uint8_t>(kShardSize, 0));
  for (std::size_t i = 0; i < kDataShardCount; ++i) {
    for (std::size_t j = 0; j < kShardSize; ++j) {
      shards[i][j] = i + j;
    }
  }
  std::cout << "Original:\n";
  Print(shards);

  if (rs.Encode(shards)) {
    std::cout << "Encoded:\n";
    Print(shards);
  }

  shards[2].clear();
  std::cout << "Clear:\n";
  Print(shards);

  if (rs.Reconstruct(shards)) {
    std::cout << "Reconstruct:\n";
    Print(shards);
  }

  std::vector<uint8_t> lost;
  std::cout << "Lost:\n";
  Print(lost);

  std::vector<std::reference_wrapper<std::vector<uint8_t>>> refs;
  refs.reserve(kTotalShardCount);
  refs.emplace_back(std::ref(lost));
  std::for_each(shards.begin() + 1, shards.end(),
                [&refs](auto& e) { refs.emplace_back(std::ref(e)); });
  if (rs.Reconstruct(refs)) {
    std::cout << "Reconstruct lost:\n";
    Print(lost);
  }

  // std::cout << "Verify: " << rs.Verify(shards) << "\n";
}
