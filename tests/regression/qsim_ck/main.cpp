#include <iostream>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <chrono>
#include <vortex.h>
#include <cmath>
#include "common.h"
#include <fstream>

#define FLOAT_ULP 6
#define MAX_NUM_QUBITS 20
#define DOUBLE_ULP 6

#define RT_CHECK(_expr)                                         \
   do {                                                         \
     int _ret = _expr;                                          \
     if (0 == _ret)                                             \
       break;                                                   \
     printf("Error: '%s' returned %d!\n", #_expr, (int)_ret);   \
	 cleanup();			                                              \
     exit(-1);                                                  \
   } while (false)

///////////////////////////////////////////////////////////////////////////////

template <typename Type>
class Comparator {};

template <>
class Comparator<int> {
public:
  static const char* type_str() {
    return "integer";
  }
  static int generate() {
    return rand();
  }
  static bool compare(int a, int b, int index, int errors) {
    if (a != b) {
      if (errors < 100) {
        printf("*** error: [%d] expected=%d, actual=%d\n", index, b, a);
      }
      return false;
    }
    return true;
  }
};

template <>
class Comparator<float> {
public:
  static const char* type_str() {
    return "float";
  }
  static int generate() {
    return static_cast<float>(rand()) / RAND_MAX;
  }
  static bool compare(float a, float b, int index, int errors) {
    union fi_t { float f; int32_t i; };
    fi_t fa, fb;
    fa.f = a;
    fb.f = b;
    auto d = std::abs(fa.i - fb.i);
    if (d > FLOAT_ULP) {
      if (errors < 100) {
        printf("*** error: [%d] expected=%f, actual=%f\n", index, b, a);
      }
      return false;
    }
    return true;
  }
};

template <>
class Comparator<complex> {
public:
  static const char* type_str() {
    return "complex";
  }
  static complex generate() {
    complex c;
    c.real = static_cast<double>(rand()) / RAND_MAX;
    c.imag = static_cast<double>(rand()) / RAND_MAX;
    return c;
  }
  static bool compare(complex a, complex b, int index, int errors) {
    union di_t { double d; int64_t i; };
    di_t dra, drb;
    di_t dia, dib;
    dra.d = a.real;
    dia.d = a.imag;
    drb.d = b.real;
    dib.d = b.imag;
    auto rd = std::abs(dra.i - drb.i);
    auto id = std::abs(dia.i - dib.i);

    if (rd > DOUBLE_ULP && id > DOUBLE_ULP) {
      if (errors < 100) {
        printf("*** error: [%d] expected=(%lf)+(%lfi), actual=(%lf)+(%lfi)\n", index, b.real, b.imag, a.real, a.imag);
      }
      return false;
    }
    return true;
  }
};

static complex mul(complex c1, complex c2) {
  complex temp;
  temp.real = c1.real * c2.real - c1.imag * c2.imag;
  temp.imag = c1.real * c2.imag - c1.imag * c2.real;
  return (temp);
}

static complex add(complex c1, complex c2) {
  complex temp;
  temp.real = c1.real + c2.real;
  temp.imag = c1.imag + c2.imag;
  return (temp);
}

static void cmatmul_cpu(TYPE* out, const TYPE* A, const TYPE* B, uint32_t width, uint32_t height) {
  for (uint32_t row = 0; row < height; ++row) {
    for (uint32_t col = 0; col < width; ++col) {
      TYPE sum = {0, 0};
      for (uint32_t e = 0; e < width; ++e) {
          sum = add(sum, mul(A[row * width + e], B[e * width + col]));
      }
      out[row * width + col] = sum;
    }
  }
}

// static void matmul_cpu(TYPE* out, const TYPE* A, const TYPE* B, uint32_t width, uint32_t height) {
//   for (uint32_t row = 0; row < height; ++row) {
//     for (uint32_t col = 0; col < width; ++col) {
//       TYPE sum(0);
//       for (uint32_t e = 0; e < width; ++e) {
//           sum += A[row * width + e] * B[e * width + col];
//       }
//       out[row * width + col] = sum;
//     }
//   }
// }

const char* kernel_file = "kernel.vxbin";
uint32_t size = 32;

vx_device_h device = nullptr;
vx_buffer_h A_buffer = nullptr;
vx_buffer_h B_buffer = nullptr;
vx_buffer_h C_buffer = nullptr;
vx_buffer_h krnl_buffer = nullptr;
vx_buffer_h args_buffer = nullptr;
kernel_arg_t kernel_arg = {};

static void show_usage() {
   std::cout << "Vortex Test." << std::endl;
   std::cout << "Usage: [-k: kernel] [-n size] [-h: help]" << std::endl;
}

static void parse_args(int argc, char **argv) {
  int c;
  while ((c = getopt(argc, argv, "n:k:h?")) != -1) {
    switch (c) {
    case 'n':
      size = atoi(optarg);
      break;
    case 'k':
      kernel_file = optarg;
      break;
    case 'h':
    case '?': {
      show_usage();
      exit(0);
    } break;
    default:
      show_usage();
      exit(-1);
    }
  }
}

void cleanup() {
  if (device) {
    vx_mem_free(A_buffer);
    vx_mem_free(B_buffer);
    vx_mem_free(C_buffer);
    vx_mem_free(krnl_buffer);
    vx_mem_free(args_buffer);
    vx_dev_close(device);
  }
}

// file -> gates vector
void parse(int num_qubits, int max_num_gates, std::string filename, std::vector<TYPE>& states, std::vector<TYPE>& gates) {
  std::ifstream infile(filename);

  // first 4 lines
  std::string temp;
  for (int i = 0; i < 8; i++) { 
    infile >> temp;
  }

  int ng_per_qubit[MAX_NUM_QUBITS] = {0, }; // #gates per qubit

  // to print gate names
  std::vector<std::vector<std::string>> gate_info(num_qubits);

  // all qubits are in 0 states
  for (int i = 0; i < num_qubits; i++) {
    complex e1 = {1, 0};
    complex e0 = {0, 0};
    states[2 * i] = e1;
    states[2 * i + 1] = e0;
  }

  std::string gate_st, qubit_st;
  while (infile >> gate_st >> qubit_st) {
      int qidx = qubit_st[qubit_st.length() - 3] - '0';
      
      // gate number
      int ng = ng_per_qubit[qidx];

      // starting index = qubit index * max #gates * 4(4 numbers in one gate)
      int startidx = qidx * max_num_gates * 4;

      gate_info[qidx].push_back(gate_st);

      if (gate_st == "h") {
        complex e1 = {0.70710678118, 0};
        complex e2 = {-0.70710678118, 0};
        gates[startidx + 4 * ng] = e1;
        gates[startidx + 4 * ng + 1] = e1;
        gates[startidx + 4 * ng + 2] = e1;
        gates[startidx + 4 * ng + 3] = e2;
      }
      else if (gate_st == "x") {
        complex e1 = {1, 0};
        complex e0 = {0, 0};
        gates[startidx + 4 * ng] = e0;
        gates[startidx + 4 * ng + 1] = e1;
        gates[startidx + 4 * ng + 2] = e1;
        gates[startidx + 4 * ng + 3] = e0;
      }
      else if (gate_st == "z") {
        complex e1 = {1, 0};
        complex e0 = {0, 0};
        complex em1 = {1, 0};
        gates[startidx + 4 * ng] = e1;
        gates[startidx + 4 * ng + 1] = e0;
        gates[startidx + 4 * ng + 2] = e0;
        gates[startidx + 4 * ng + 3] = em1;
      }
    
    ng_per_qubit[qidx]++;
  }

  // store 100(flag) at the end of each gate array
  for (uint32_t i = 0; i < num_qubits; i++) {
    int startidx = 4 * max_num_gates * i;
    int ng = ng_per_qubit[i];
    complex t = {100, 0};
    gates[startidx + 4 * ng] = t;    // end flag
  }

  std::cout << "\n\nAfter Parsing" << "\n";
  for (uint32_t i = 0; i < num_qubits; i++) {
    std::cout << "Qubit-" << i << "  ";
    for (uint32_t g = 0; g < ng_per_qubit[i]; g++) {
      std::cout << gate_info[i][g] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "\n\n";
}

int main(int argc, char *argv[]) {
  // parse command arguments
  parse_args(argc, argv);

  // open device connection
  std::cout << "open device connection" << std::endl;
  RT_CHECK(vx_dev_open(&device));

  // std::srand(50);
  std::string filename = "/home/ubuntu/testq/qsim_vortex/tests/regression/qsim2/temp.qasm";
  std::ifstream infile(filename);

  // first 4 lines
  std::string temp;
  for (int i = 0; i < 8; i++) { 
    infile >> temp;
  }
  uint32_t num_qubits = temp[temp.length() - 3] - '0'; 
  uint32_t max_num_gates = 20;


  // generate source data
  // store state vector of each qubit
  std::vector<TYPE> h_A(2 * num_qubits);
  // store gate array on each qubit
  std::vector<TYPE> h_B(4 * max_num_gates * num_qubits);
  // result state vector
  std::vector<TYPE> h_C(2 * num_qubits);
  
  parse(num_qubits, max_num_gates, filename, h_A, h_B);

  uint32_t states_buf_size = 2 * num_qubits * sizeof(TYPE);
  uint32_t gates_buf_size = 4 * max_num_gates * num_qubits * sizeof(TYPE);

  std::cout << "data type: " << Comparator<TYPE>::type_str() << std::endl;
  std::cout << "nq:" << num_qubits << " sbs: " << states_buf_size << " gbs: " << gates_buf_size << std::endl;

  kernel_arg.num_tasks = num_qubits;
  kernel_arg.max_num_gates = max_num_gates;
  kernel_arg.nq = num_qubits;

  // allocate device memory
  std::cout << "allocate device memory" << std::endl;
  RT_CHECK(vx_mem_alloc(device, states_buf_size, VX_MEM_READ, &A_buffer));
  RT_CHECK(vx_mem_address(A_buffer, &kernel_arg.A_addr));
  RT_CHECK(vx_mem_alloc(device, gates_buf_size, VX_MEM_READ, &B_buffer));
  RT_CHECK(vx_mem_address(B_buffer, &kernel_arg.B_addr));
  RT_CHECK(vx_mem_alloc(device, states_buf_size, VX_MEM_WRITE, &C_buffer));
  RT_CHECK(vx_mem_address(C_buffer, &kernel_arg.C_addr));

  std::cout << "A_addr=0x" << std::hex << kernel_arg.A_addr << std::endl;
  std::cout << "B_addr=0x" << std::hex << kernel_arg.B_addr << std::endl;
  std::cout << "C_addr=0x" << std::hex << kernel_arg.C_addr << std::endl;

  // upload matrix A buffer
  {
    std::cout << "upload matrix A buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(A_buffer, h_A.data(), 0, states_buf_size));
  }

  // upload matrix B buffer
  {
    std::cout << "upload matrix B buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(B_buffer, h_B.data(), 0, gates_buf_size));
  }

  // upload program
  std::cout << "upload program" << std::endl;
  RT_CHECK(vx_upload_kernel_file(device, kernel_file, &krnl_buffer));

  // upload kernel argument
  std::cout << "upload kernel argument" << std::endl;
  RT_CHECK(vx_upload_bytes(device, &kernel_arg, sizeof(kernel_arg_t), &args_buffer));

  auto time_start = std::chrono::high_resolution_clock::now();

  // start device
  std::cout << "start device" << std::endl;
  RT_CHECK(vx_start(device, krnl_buffer, args_buffer));

  // wait for completion
  std::cout << "wait for completion" << std::endl;
  RT_CHECK(vx_ready_wait(device, VX_MAX_TIMEOUT));

  auto time_end = std::chrono::high_resolution_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
  printf("Elapsed time: %lg ms\n", elapsed);

  // download destination buffer
  std::cout << "download destination buffer" << std::endl;
  RT_CHECK(vx_copy_from_dev(h_C.data(), C_buffer, 0, states_buf_size));


  std::cout << "\n\n";
  for (uint32_t i = 0; i < num_qubits; i++) {
    std::cout << "qubit " << i << " ";
    std::cout << "(" << h_C[2 * i].real << ")+(" << h_C[2 * i].imag << "i) ";
    std::cout << "(" << h_C[2 * i + 1].real << ")+(" << h_C[2 * i + 1].imag << "i) " << std::endl;
  }
  std::cout << "\n\n";


  // cleanup
  std::cout << "cleanup" << std::endl;
  cleanup();

  std::cout << "PASSED!" << std::endl;

  return 0;
}