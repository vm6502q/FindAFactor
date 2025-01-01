//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "oclengine.hpp"

#include <algorithm>
#include <iostream>
#include <regex>
#include <sstream>

namespace Qimcifa {

/// "Qrack::OCLEngine" manages the single OpenCL context

// Public singleton methods to get pointers to various methods
DeviceContextPtr OCLEngine::GetDeviceContextPtr(const int64_t& dev)
{
    if ((dev >= GetDeviceCount()) || (dev < -1) || (dev >= ((int64_t)all_device_contexts.size()))) {
        throw std::invalid_argument("Invalid OpenCL device selection");
    } else if (dev == -1) {
        return default_device_context;
    } else {
        return all_device_contexts[dev];
    }
}

// clang-format off
const std::vector<OCLKernelHandle> OCLEngine::kernelHandles{
    OCLKernelHandle(OCL_API_FACTORIZE_SMOOTH, "factorize")
};
// clang-format on

const std::string OCLEngine::binary_file_prefix("ocl_dev_");
const std::string OCLEngine::binary_file_ext(".ir");

std::vector<DeviceContextPtr> OCLEngine::GetDeviceContextPtrVector() { return all_device_contexts; }
void OCLEngine::SetDeviceContextPtrVector(std::vector<DeviceContextPtr> vec, DeviceContextPtr dcp)
{
    all_device_contexts = vec;
    if (dcp != nullptr) {
        default_device_context = dcp;
    }
}

void OCLEngine::SetDefaultDeviceContext(DeviceContextPtr dcp) { default_device_context = dcp; }

cl::Program OCLEngine::MakeProgram(std::shared_ptr<OCLDeviceContext> devCntxt)
{
    // Load and build kernel
    std::string kernelSourceStr =
        "#define BCAPPOW " + std::string(getenv("BCAPPOW")) + "\n" +
        "#define BIG_INTEGER_WORD_BITS 64U\n" +
        "#define BIG_INTEGER_WORD_POWER 6U\n" +
        "#define BIG_INTEGER_WORD uint64_t\n" +
        "#define BIG_INTEGER_HALF_WORD uint32_t\n" +
        "#define BIG_INTEGER_HALF_WORD_MASK 0xFFFFFFFFULL\n" +
        "#define BIG_INTEGER_HALF_WORD_MASK_NOT 0xFFFFFFFF00000000ULL\n" +
        "\n" +
        "// This can be any power of 2 greater than (or equal to) 64:\n" +
        "const size_t BIG_INTEGER_BITS = (1 << BCAPPOW);\n" +
        "const int BIG_INTEGER_WORD_SIZE = BIG_INTEGER_BITS / BIG_INTEGER_WORD_BITS;\n" +
        "\n" +
        "// The rest of the constants need to be consistent with the one above:\n" +
        "const size_t BIG_INTEGER_HALF_WORD_BITS = BIG_INTEGER_WORD_BITS >> 1U;\n" +
        "const int BIG_INTEGER_HALF_WORD_SIZE = BIG_INTEGER_WORD_SIZE << 1U;\n" +
        "const int BIG_INTEGER_MAX_WORD_INDEX = BIG_INTEGER_WORD_SIZE - 1U;\n" +
        "\n" +
        "typedef struct BigInteger {\n" +
        "    BIG_INTEGER_WORD bits[BIG_INTEGER_WORD_SIZE];\n" +
        "\n" +
        "    inline BigInteger()\n" +
        "    {\n" +
        "        // Intentionally left blank.\n" +
        "    }\n" +
        "\n" +
        "    inline BigInteger(const BigInteger& val)\n" +
        "    {\n" +
        "        for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {\n" +
        "            bits[i] = val.bits[i];\n" +
        "        }\n" +
        "    }\n" +
        "\n" +
        "    inline BigInteger(const BIG_INTEGER_WORD& val)\n" +
        "    {\n" +
        "        bits[0] = val;\n" +
        "        for (int i = 1; i < BIG_INTEGER_WORD_SIZE; ++i) {\n" +
        "            bits[i] = 0U;\n" +
        "        }\n" +
        "    }\n" +
        "\n" +
        "    inline void set_0()\n" +
        "    {\n" +
        "        for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {\n" +
        "            bits[i] = 0U;\n" +
        "        }\n" +
        "    }\n" +
        "\n" +
        "    inline void xor_bit(const BIG_INTEGER_HALF_WORD& b) {\n" +
        "        bits[b % BIG_INTEGER_WORD_BITS] ^= (1ULL << (b / BIG_INTEGER_WORD_BITS));\n" +
        "    }\n" +
        "} BigInteger;\n" +
        "\n" +
        "// \"Schoolbook division\" (on half words)\n" +
        "// Complexity - O(x^2)\n" +
        "void bi_div_mod_small(\n" +
        "    const BigInteger& left, BIG_INTEGER_HALF_WORD right, BigInteger* quotient, BIG_INTEGER_HALF_WORD* rmndr)\n" +
        "{\n" +
        "    BIG_INTEGER_WORD carry = 0U;\n" +
        "    if (quotient) {\n" +
        "        quotient.set_0();\n" +
        "        for (int i = BIG_INTEGER_HALF_WORD_SIZE - 1; i >= 0; --i) {\n" +
        "            const int i2 = i >> 1;\n" +
        "            carry <<= BIG_INTEGER_HALF_WORD_BITS;\n" +
        "            if (i & 1) {\n" +
        "                carry |= left.bits[i2] >> BIG_INTEGER_HALF_WORD_BITS;\n" +
        "                quotient->bits[i2] |= (carry / right) << BIG_INTEGER_HALF_WORD_BITS;\n" +
        "            } else {\n" +
        "                carry |= left.bits[i2] & BIG_INTEGER_HALF_WORD_MASK;\n" +
        "                quotient->bits[i2] |= (carry / right);\n" +
        "            }\n" +
        "            carry %= right;\n" +
        "        }\n" +
        "    } else {\n" +
        "        for (int i = BIG_INTEGER_HALF_WORD_SIZE - 1; i >= 0; --i) {\n" +
        "            const int i2 = i >> 1;\n" +
        "            carry <<= BIG_INTEGER_HALF_WORD_BITS;\n" +
        "            if (i & 1) {\n" +
        "                carry |= left.bits[i2] >> BIG_INTEGER_HALF_WORD_BITS;\n" +
        "            } else {\n" +
        "                carry |= left.bits[i2] & BIG_INTEGER_HALF_WORD_MASK;\n" +
        "            }\n" +
        "            carry %= right;\n" +
        "        }\n" +
        "    }\n" +
        "\n" +
        "    *rmndr = carry;\n" +
        "}\n" +
        "\n" +
        "__kernel void factorize(\n" +
        "    __global const BigInteger *numbers,                    // Array of numbers to check\n" +
        "    __global const int *primes,                            // Array of small primes for smoothness\n" +
        "    __global bool *results,                                // Output: 1 if smooth, 0 if not\n" +
        "    __global const BigInteger *factor_vectors,             // Output: Factorization vectors as bitmasks\n" +
        "    const int primeCount                                   // Number of primes in the array\n" +
        ") {\n" +
        "    int gid = get_global_id(0);                            // Get the index of this work item\n" +
        "    BigInteger number = numbers[gid];                      // The number to check\n" +
        "    BigInteger factor_vector = 0U;                         // Initialize the factor vector as 0\n" +
        "    BigInteger q;                                          // For quotient\n" +
        "\n" +
        "    // Test divisibility by each prime\n" +
        "    for (int i = 0; i < primeCount; ++i) {\n" +
        "        const int& p = primes[i];\n" +
        "        do {\n" +
        "            unsigned int r = 0U;\n" +
        "            bi_div_mod_small(number, primes[i], &q, &r);\n" +
        "            if (!r) {\n" +
        "                number = q;\n" +
        "                factor_vector.xor_bit(i);                  // Flip the corresponding bit\n" +
        "            }\n" +
        "        } while (number >= p)\n" +
        "    }\n" +
        "\n" +
        "    // If number is reduced to 1, it is smooth\n" +
        "    results[gid] = (number == 1);\n" +
        "\n" +
        "    // Store the factor vector\n" +
        "    factor_vectors[gid] = factor_vector;\n" +
        "}\n";

    cl::Program::Sources sources;
    sources.push_back({(const char*)kernelSourceStr.c_str(), (long unsigned int)(kernelSourceStr.size() + 1U) });

    cl::Program program = cl::Program(devCntxt->context, sources);
    std::cout << "Building JIT." << std::endl;

    return program;
}

void OCLEngine::SaveBinary(cl::Program program, std::string path, std::string fileName)
{
    std::vector<size_t> clBinSizes = program.getInfo<CL_PROGRAM_BINARY_SIZES>();
    size_t clBinSize = 0U;
    int64_t clBinIndex = 0;

    for (size_t i = 0U; i < clBinSizes.size(); ++i) {
        if (clBinSizes[i]) {
            clBinSize = clBinSizes[i];
            clBinIndex = i;
            break;
        }
    }

    std::cout << "Binary size:" << clBinSize << std::endl;

#if defined(_WIN32) && !defined(__CYGWIN__)
    int err = _mkdir(path.c_str());
#else
    int err = mkdir(path.c_str(), 0700);
#endif
    if (err != -1) {
        std::cout << "Making directory: " << path << std::endl;
    }

    FILE* clBinFile = fopen((path + fileName).c_str(), "w");
    std::vector<std::vector<unsigned char>> clBinaries = program.getInfo<CL_PROGRAM_BINARIES>();
    std::vector<unsigned char> clBinary = clBinaries[clBinIndex];
    fwrite(&clBinary[0U], clBinSize, sizeof(unsigned char), clBinFile);
    fclose(clBinFile);
}

InitOClResult OCLEngine::InitOCL(std::vector<int64_t> maxAllocVec)
{
    // get all platforms (drivers), e.g. NVIDIA

    std::vector<cl::Platform> all_platforms;
    std::vector<cl::Device> all_devices;
    std::vector<int64_t> device_platform_id;
    cl::Platform default_platform;
    cl::Device default_device;
    std::vector<DeviceContextPtr> all_dev_contexts;
    DeviceContextPtr default_dev_context;

    cl::Platform::get(&all_platforms);

    if (all_platforms.empty()) {
        std::cout << " No platforms found. Check OpenCL installation!\n";
        return InitOClResult();
    }

    // get all devices
    std::vector<cl::Platform> devPlatVec;
    std::vector<std::vector<cl::Device>> all_platforms_devices;
    std::vector<bool> all_devices_is_gpu;
    std::vector<bool> all_devices_is_cpu;
    for (size_t i = 0U; i < all_platforms.size(); ++i) {
        all_platforms_devices.push_back(std::vector<cl::Device>());
        all_platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &(all_platforms_devices[i]));
        for (size_t j = 0U; j < all_platforms_devices[i].size(); ++j) {
            // VirtualCL seems to break if the assignment constructor of cl::Platform is used here from the original
            // list. Assigning the object from a new query is always fine, though. (They carry the same underlying
            // platform IDs.)
            std::vector<cl::Platform> temp_platforms;
            cl::Platform::get(&temp_platforms);
            devPlatVec.push_back(temp_platforms[i]);
            device_platform_id.push_back(i);
        }
        all_devices.insert(all_devices.end(), all_platforms_devices[i].begin(), all_platforms_devices[i].end());

        // Linux implements `cl::Device` relation operators, including equality, but Mac considers OpenCL "deprecated,"
        // and other compilers might not see a strict need in OpenCL implementation standard for a `cl::Device` equality
        // operator, which would allow the use of `std::find()`.
        std::vector<cl::Device> gpu_devices;
        all_platforms[i].getDevices(CL_DEVICE_TYPE_GPU, &gpu_devices);
        std::vector<bool> gpu_to_insert(all_platforms_devices[i].size(), false);
        for (size_t j = 0U; j < gpu_devices.size(); ++j) {
            for (size_t k = 0U; k < all_platforms_devices[i].size(); ++k) {
                if (gpu_devices[j].getInfo<CL_DEVICE_NAME>() == all_platforms_devices[i][j].getInfo<CL_DEVICE_NAME>()) {
                    // Assuming all devices with the same name are identical vendor, line, and model, this works.
                    gpu_to_insert[k] = true;
                }
            }
        }
        all_devices_is_gpu.insert(all_devices_is_gpu.end(), gpu_to_insert.begin(), gpu_to_insert.end());

        std::vector<cl::Device> cpu_devices;
        all_platforms[i].getDevices(CL_DEVICE_TYPE_CPU, &cpu_devices);
        std::vector<bool> cpu_to_insert(all_platforms_devices[i].size(), false);
        for (size_t j = 0U; j < cpu_devices.size(); ++j) {
            for (size_t k = 0U; k < all_platforms_devices[i].size(); ++k) {
                if (cpu_devices[j].getInfo<CL_DEVICE_NAME>() == all_platforms_devices[i][j].getInfo<CL_DEVICE_NAME>()) {
                    // Assuming all devices with the same name are identical vendor, line, and model, this works.
                    cpu_to_insert[k] = true;
                }
            }
        }
        all_devices_is_cpu.insert(all_devices_is_cpu.end(), cpu_to_insert.begin(), cpu_to_insert.end());
    }
    if (all_devices.empty()) {
        std::cout << " No devices found. Check OpenCL installation!\n";
        return InitOClResult();
    }

    int64_t deviceCount = all_devices.size();
    // prefer the last device because that's usually a GPU or accelerator; device[0U] is usually the CPU
    int64_t dev = deviceCount - 1;
    if (getenv("FINDAFACTOR_OCL_DEFAULT_DEVICE")) {
        dev = std::stoi(std::string(getenv("FINDAFACTOR_OCL_DEFAULT_DEVICE")));
        if ((dev < 0) || (dev > (deviceCount - 1))) {
            std::cout << "WARNING: Invalid FINDAFACTOR_OCL_DEFAULT_DEVICE selection. (Falling back to highest index device "
                         "as default.)"
                      << std::endl;
            dev = deviceCount - 1;
        }
    }

    // create the programs that we want to execute on the devices
    int64_t plat_id = -1;
    std::vector<cl::Context> all_contexts;
    std::vector<std::string> all_filenames;
    for (int64_t i = 0; i < deviceCount; ++i) {
        // a context is like a "runtime link" to the device and platform;
        // i.e. communication is possible
        if (device_platform_id[i] != plat_id) {
            plat_id = device_platform_id[i];
            all_contexts.push_back(cl::Context(all_platforms_devices[plat_id]));
        }
        const std::string devName(all_devices[i].getInfo<CL_DEVICE_NAME>());
        const bool useHostRam = all_devices_is_cpu[i] || (devName.find("Intel(R) UHD") != std::string::npos) ||
            (devName.find("Iris") != std::string::npos);
        DeviceContextPtr devCntxt =
            std::make_shared<OCLDeviceContext>(devPlatVec[i], all_devices[i], all_contexts[all_contexts.size() - 1U], i,
                plat_id, maxAllocVec[i % maxAllocVec.size()], all_devices_is_gpu[i], all_devices_is_cpu[i], useHostRam);

        std::cout << "Device #" << i << ", ";
        cl::Program program = MakeProgram(devCntxt);

        cl_int buildError =
            program.build({ all_devices[i] }, "-cl-strict-aliasing -cl-denorms-are-zero -cl-fast-relaxed-math");
        if (buildError != CL_SUCCESS) {
            std::cout << "Error building for device #" << i << ": " << buildError << ", "
                      << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(all_devices[i])
                      << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(all_devices[i]) << std::endl;

            // The default device was set above to be the last device in the list. If we can't compile for it, we
            // use the first device. If the default is the first device, and we can't compile for it, then we don't
            // have any devices that can compile at all, and the environment needs to be fixed by the user.
            if (i == dev) {
                default_dev_context = all_dev_contexts[0U];
                default_platform = all_platforms[0U];
                default_device = all_devices[0U];
            }

            continue;
        }

        all_dev_contexts.push_back(devCntxt);

        for (unsigned int j = 0U; j < kernelHandles.size(); ++j) {
            all_dev_contexts[i]->calls[kernelHandles[j].oclapi] =
                cl::Kernel(program, kernelHandles[j].kernelname.c_str());
            all_dev_contexts[i]->mutexes.emplace(kernelHandles[j].oclapi, new std::mutex);
        }

        if (i == dev) {
            default_dev_context = all_dev_contexts[i];
            default_platform = all_platforms[plat_id];
            default_device = all_devices[i];
        }
    }

    // For VirtualCL support, the device info can only be accessed AFTER all contexts are created.
    std::cout << "Default platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";
    std::cout << "Default device: #" << dev << ", " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";
    for (int64_t i = 0; i < deviceCount; ++i) {
        std::cout << "OpenCL device #" << i << ": " << all_devices[i].getInfo<CL_DEVICE_NAME>() << "\n";
    }

    return InitOClResult(all_dev_contexts, default_dev_context);
}

OCLEngine::OCLEngine()
    : maxActiveAllocSizes(1U, -1)
{
    if (getenv("FINDAFACTOR_MAX_ALLOC_MB")) {
        std::string devListStr = std::string(getenv("FINDAFACTOR_MAX_ALLOC_MB"));
        maxActiveAllocSizes.clear();
        if (devListStr.compare("")) {
            std::stringstream devListStr_stream(devListStr);
            // See
            // https://stackoverflow.com/questions/7621727/split-a-string-into-words-by-multiple-delimiters#answer-58164098
            std::regex re("[.]");
            while (devListStr_stream.good()) {
                std::string term;
                getline(devListStr_stream, term, ',');
                // the '-1' is what makes the regex split (-1 := what was not matched)
                std::sregex_token_iterator first{ term.begin(), term.end(), re, -1 }, last;
                std::vector<std::string> tokens{ first, last };
                if (tokens.size() == 1U) {
                    maxActiveAllocSizes.push_back(stoi(term));
                    if (maxActiveAllocSizes.back() >= 0) {
                        maxActiveAllocSizes.back() = maxActiveAllocSizes.back() << 20U;
                    }
                    continue;
                }
                const unsigned maxI = stoi(tokens[0U]);
                std::vector<int64_t> limits(tokens.size() - 1U);
                for (unsigned i = 1U; i < tokens.size(); ++i) {
                    limits[i - 1U] = stoi(tokens[i]);
                }
                for (unsigned i = 0U; i < maxI; ++i) {
                    for (unsigned j = 0U; j < limits.size(); ++j) {
                        maxActiveAllocSizes.push_back(limits[j]);
                        if (maxActiveAllocSizes.back() >= 0) {
                            maxActiveAllocSizes.back() = maxActiveAllocSizes.back() << 20U;
                        }
                    }
                }
            }
        }
    }

    InitOClResult initResult = InitOCL(maxActiveAllocSizes);
    SetDeviceContextPtrVector(initResult.all_dev_contexts, initResult.default_dev_context);
    activeAllocSizes = std::vector<size_t>(initResult.all_dev_contexts.size());
}

} // namespace Qrack
