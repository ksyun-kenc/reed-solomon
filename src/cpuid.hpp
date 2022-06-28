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

#if __GNUC__
#include <cpuid.h>
#elif _MSC_VER
#include <intrin.h>
#else
#error Only supports CLANG/GCC or MSVC.
#endif

#include <array>
#include <bitset>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// eax_value: info type/function ID
void GetCpuIdEx(int32_t cpu_info[4], int32_t eax_value, int32_t ecx_value) {
#if __GNUC__
  __cpuid_count(eax_value, ecx_value, cpu_info[0], cpu_info[1], cpu_info[2],
                cpu_info[3]);
#elif _MSC_VER
#if _WIN64 || _MSC_VER >= 1600
  __cpuidex(cpu_info, eax_value, ecx_value);
#else
  if (nullptr == cpu_info) {
    return;
  }
  _asm {
    mov edi, cpu_info;
    mov eax, eax_value;
    mov ecx, ecx_value;
    cpuid;
    mov    [edi], eax;
    mov    [edi+4], ebx;
    mov    [edi+8], ecx;
    mov    [edi+12], edx;
  }
#endif
#endif
}

void GetCpuId(int32_t cpu_info[4], int32_t eax_value) {
#if __GNUC__
  __cpuid(eax_value, cpu_info[0], cpu_info[1], cpu_info[2], cpu_info[3]);
#elif _MSC_VER
#if _MSC_VER >= 1400
  __cpuid(cpu_info, eax_value);
#else
  GetCpuIdEx(cpu_info, eax_value, 0);
#endif
#endif
}

class CpuId {
  class CpuIdInternal;

 public:
  static std::string Vendor(void) { return cpu_info_.vendor_; }
  static std::string Brand(void) { return cpu_info_.brand_; }

  static bool SSE3(void) { return cpu_info_.f_1_0_ecx_[0]; }
  static bool PCLMULQDQ(void) { return cpu_info_.f_1_0_ecx_[1]; }
  static bool MONITOR(void) { return cpu_info_.f_1_0_ecx_[3]; }
  static bool SSSE3(void) { return cpu_info_.f_1_0_ecx_[9]; }
  static bool FMA(void) { return cpu_info_.f_1_0_ecx_[12]; }
  static bool FMA3(void) { return FMA() && OSXSAVE(); }
  static bool CMPXCHG16B(void) { return cpu_info_.f_1_0_ecx_[13]; }
  static bool SSE41(void) { return cpu_info_.f_1_0_ecx_[19]; }
  static bool SSE42(void) { return cpu_info_.f_1_0_ecx_[20]; }
  static bool MOVBE(void) { return cpu_info_.f_1_0_ecx_[22]; }
  static bool POPCNT(void) { return cpu_info_.f_1_0_ecx_[23]; }
  static bool AES(void) { return cpu_info_.f_1_0_ecx_[25]; }
  static bool XSAVE(void) { return cpu_info_.f_1_0_ecx_[26]; }
  static bool OSXSAVE(void) { return cpu_info_.f_1_0_ecx_[27]; }
  static bool AVX(void) { return cpu_info_.f_1_0_ecx_[28]; }
  static bool F16C(void) { return cpu_info_.f_1_0_ecx_[29]; }
  static bool RDRAND(void) { return cpu_info_.f_1_0_ecx_[30]; }

  static bool MSR(void) { return cpu_info_.f_1_0_edx_[5]; }
  static bool CX8(void) { return cpu_info_.f_1_0_edx_[8]; }
  static bool SEP(void) { return cpu_info_.f_1_0_edx_[11]; }
  static bool CMOV(void) { return cpu_info_.f_1_0_edx_[15]; }
  static bool CLFSH(void) { return cpu_info_.f_1_0_edx_[19]; }
  static bool MMX(void) { return cpu_info_.f_1_0_edx_[23]; }
  static bool FXSR(void) { return cpu_info_.f_1_0_edx_[24]; }
  static bool SSE(void) { return cpu_info_.f_1_0_edx_[25]; }
  static bool SSE2(void) { return cpu_info_.f_1_0_edx_[26]; }

  static bool FSGSBASE(void) { return cpu_info_.f_7_0_ebx_[0]; }
  static bool BMI1(void) { return cpu_info_.f_7_0_ebx_[3]; }
  static bool HLE(void) {
    return cpu_info_.is_intel_ && cpu_info_.f_7_0_ebx_[4];
  }
  static bool AVX2(void) { return cpu_info_.f_7_0_ebx_[5]; }
  static bool BMI2(void) { return cpu_info_.f_7_0_ebx_[8]; }
  static bool ERMS(void) { return cpu_info_.f_7_0_ebx_[9]; }
  static bool INVPCID(void) { return cpu_info_.f_7_0_ebx_[10]; }
  static bool RTM(void) {
    return cpu_info_.is_intel_ && cpu_info_.f_7_0_ebx_[11];
  }
  static bool AVX512F(void) { return cpu_info_.f_7_0_ebx_[16]; }
  static bool AVX512DQ(void) { return cpu_info_.f_7_0_ebx_[17]; }
  static bool RDSEED(void) { return cpu_info_.f_7_0_ebx_[18]; }
  static bool ADX(void) { return cpu_info_.f_7_0_ebx_[19]; }
  static bool AVX512IFMA(void) { return cpu_info_.f_7_0_ebx_[21]; }
  static bool AVX512PF(void) { return cpu_info_.f_7_0_ebx_[26]; }
  static bool AVX512ER(void) { return cpu_info_.f_7_0_ebx_[27]; }
  static bool AVX512CD(void) { return cpu_info_.f_7_0_ebx_[28]; }
  static bool SHA(void) { return cpu_info_.f_7_0_ebx_[29]; }
  static bool AVX512BW(void) { return cpu_info_.f_7_0_ebx_[30]; }
  static bool AVX512VL(void) { return cpu_info_.f_7_0_ebx_[31]; }

  static bool PREFETCHWT1(void) { return cpu_info_.f_7_0_ecx_[0]; }
  static bool AVX512VBMI(void) { return cpu_info_.f_7_0_ecx_[1]; }
  static bool AVX512VBMI2(void) { return cpu_info_.f_7_0_ecx_[6]; }
  static bool GFNI(void) { return cpu_info_.f_7_0_ecx_[8]; }
  static bool VAES(void) { return cpu_info_.f_7_0_ecx_[9]; }
  static bool VPCLMULQDQ(void) { return cpu_info_.f_7_0_ecx_[10]; }
  static bool AVX512VNNI(void) { return cpu_info_.f_7_0_ecx_[11]; }
  static bool AVX512BITALG(void) { return cpu_info_.f_7_0_ecx_[12]; }
  static bool AVX512VPOPCNTDQ(void) { return cpu_info_.f_7_0_ecx_[14]; }

  static bool AVX512VP2INTERSECT(void) { return cpu_info_.f_7_0_edx_[8]; }
  static bool AMXBF16(void) { return cpu_info_.f_7_0_edx_[22]; }
  static bool AMXTILE(void) { return cpu_info_.f_7_0_edx_[24]; }
  static bool AMXINT8(void) { return cpu_info_.f_7_0_edx_[25]; }

  static bool AVX512BF16(void) { return cpu_info_.f_7_1_eax_[5]; }

  static bool LAHF(void) { return cpu_info_.f_81_0_ecx_[0]; }
  static bool LZCNT(void) {
    return cpu_info_.is_intel_ && cpu_info_.f_81_0_ecx_[5];
  }
  static bool ABM(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_ecx_[5];
  }
  static bool SSE4a(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_ecx_[6];
  }
  static bool XOP(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_ecx_[11];
  }
  static bool FMA4(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_ecx_[17];
  }
  static bool TBM(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_ecx_[21];
  }

  static bool SYSCALL(void) {
    return cpu_info_.is_intel_ && cpu_info_.f_81_0_edx_[11];
  }
  static bool MMXEXT(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_edx_[22];
  }
  static bool RDTSCP(void) {
    return cpu_info_.is_intel_ && cpu_info_.f_81_0_edx_[27];
  }
  static bool _3DNOWEXT(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_edx_[30];
  }
  static bool _3DNOW(void) {
    return cpu_info_.is_amd_ && cpu_info_.f_81_0_edx_[31];
  }

 private:
  static const CpuIdInternal cpu_info_;

  class CpuIdInternal {
   public:
    CpuIdInternal()
        : highest_id_{0},
          highest_ext_id_{0},
          is_intel_{false},
          is_amd_{false},
          f_1_0_ecx_{0},
          f_1_0_edx_{0},
          f_7_0_ebx_{0},
          f_7_0_ecx_{0},
          f_7_0_edx_{0},
          f_7_1_eax_{0},
          f_81_0_ecx_{0},
          f_81_0_edx_{0} {
      std::array<int, 4> cpui;

      // Calling GetCpuId with 0x0 as the function_id argument gets highest
      // available non-extended function_id value supported by the processor.
      GetCpuId(cpui.data(), 0);
      highest_id_ = cpui[0];

      // Capture vendor string
      char vendor[0x20] = {};
      GetCpuIdEx(cpui.data(), 0, 0);
      *reinterpret_cast<int*>(vendor) = cpui[1];
      *reinterpret_cast<int*>(vendor + 4) = cpui[3];
      *reinterpret_cast<int*>(vendor + 8) = cpui[2];
      vendor_ = vendor;
      if (vendor_ == "GenuineIntel") {
        is_intel_ = true;
      } else if (vendor_ == "AuthenticAMD") {
        is_amd_ = true;
      }

      // load bitset with flags for function 0x00000001
      if (1 <= highest_id_) {
        GetCpuIdEx(cpui.data(), 1, 0);
        f_1_0_ecx_ = cpui[2];
        f_1_0_edx_ = cpui[3];
      }

      // load bitset with flags for function 0x00000007
      if (7 <= highest_id_) {
        GetCpuIdEx(cpui.data(), 7, 0);
        f_7_0_ebx_ = cpui[1];
        f_7_0_ecx_ = cpui[2];
        f_7_0_edx_ = cpui[3];

        GetCpuIdEx(cpui.data(), 7, 1);
        f_7_1_eax_ = cpui[0];
      }

      // Calling GetCpuId with 0x80000000 as the function_id argument gets the
      // highest valid extended ID.
      GetCpuId(cpui.data(), 0x80000000);
      highest_ext_id_ = cpui[0];

      // load bitset with flags for function 0x80000001
      if (0x80000001 <= highest_ext_id_) {
        GetCpuIdEx(cpui.data(), 0x80000001, 0);
        f_81_0_ecx_ = cpui[2];
        f_81_0_edx_ = cpui[3];
      }

      // Interpret CPU brand string if reported
      char brand[0x40] = {};
      if (0x80000004 <= highest_ext_id_) {
        GetCpuIdEx(cpui.data(), 0x80000002, 0);
        std::memcpy(brand, cpui.data(), sizeof(cpui));
        GetCpuIdEx(cpui.data(), 0x80000003, 0);
        std::memcpy(brand + 16, cpui.data(), sizeof(cpui));
        GetCpuIdEx(cpui.data(), 0x80000004, 0);
        std::memcpy(brand + 32, cpui.data(), sizeof(cpui));
        brand_ = brand;
      }
    };

    unsigned int highest_id_;
    unsigned int highest_ext_id_;
    std::string vendor_;
    std::string brand_;
    bool is_intel_;
    bool is_amd_;

    std::bitset<32> f_1_0_ecx_;
    std::bitset<32> f_1_0_edx_;

    std::bitset<32> f_7_0_ebx_;
    std::bitset<32> f_7_0_ecx_;
    std::bitset<32> f_7_0_edx_;

    std::bitset<32> f_7_1_eax_;

    std::bitset<32> f_81_0_ecx_;
    std::bitset<32> f_81_0_edx_;
  };
};