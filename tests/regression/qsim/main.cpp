#include <iostream>
#include <unistd.h>
#include <string.h>
#include <vector>
#include <chrono>
#include <vortex.h>
#include <cmath>
#include "common.h"
#include "gates.h"
#include <fstream>
#include <sstream>

#define FLOAT_ULP 6
#define DOUBLE_ULP 6
#define MAX_NUM_QUBITS 20

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
class Comparator<double> {
public:
  static const char* type_str() {
    return "double";
  }
  // static double generate() {
  //   double c;
  //   c.real = static_cast<double>(rand()) / RAND_MAX;
  //   c.imag = static_cast<double>(rand()) / RAND_MAX;
  //   return c;
  // }
  static bool compare(double a, double b, int index, int errors) {
    idouble da, db;
    ifloat ra, rb, ia, ib;

    da.d = a;
    db.d = b;

    ra.u32 = (da.u64 >> 32);
    ia.u32 = (da.u64 & 0xffffffff);

    rb.u32 = (db.u64 >> 32);
    ib.u32 = (db.u64 & 0xffffffff);

    auto rd = std::abs(ra.i32 - rb.i32);
    auto id = std::abs(ia.i32 - ib.i32);

    if (rd > DOUBLE_ULP && id > DOUBLE_ULP) {
      if (errors < 100) {
        printf("*** error: [%d] expected=(%f)+(%fi), actual=(%f)+(%fi)\n", index, ra.f, ia.f, rb.f, ib.f);
      }
      return false;
    }
    return true;
  }
};

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
vx_buffer_h state_buffer = nullptr;
vx_buffer_h nqubit_buffer = nullptr;
vx_buffer_h ind_buffer = nullptr;
vx_buffer_h matrix_buffer = nullptr;
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
    vx_mem_free(state_buffer);
    vx_mem_free(nqubit_buffer);
    vx_mem_free(ind_buffer);
    vx_mem_free(matrix_buffer);
    vx_mem_free(krnl_buffer);
    vx_mem_free(args_buffer);
    vx_dev_close(device);
  }
}

// void generate_gates(std::string gate, std::vector<int>& qubits, std::vector<Gate>& gates);

Gate* generate_gates(std::string gate, std::vector<int>& qubits, std::vector<Gate*>& gates) {
    Gate* g = nullptr;
    if (gate == "h")
      g = new HGate();
    else if (gate == "x")
      g = new XGate();
    else if (gate == "z")
      g = new ZGate();
    else if (gate == "t")
      g = new TGate();
    else if (gate == "tdg")
      g = new TdgGate();
    else if (gate == "cx")
      g = new CxGate();

    if (g != nullptr) {
      g->initialize(qubits);
      // gates.push_back(g);
    }
    return g;
}

void convert_gates_to_data(std::vector<Gate*>& gates, std::vector<int>& num_qubits, 
  std::vector<int>& qubit_indexes, std::vector<TYPE>& gate_matrices) {
  uint32_t q_ind = 0;
  uint32_t ind_ind = 0;
  uint32_t mat_ind = 0;
  for(Gate* g : gates) {
    num_qubits[q_ind++] = (g->qubits.size());
    
    for(int qubit_index : g->qubits)
      qubit_indexes[ind_ind++] = (qubit_index);
    
    for(TYPE data: g->matrix)
      gate_matrices[mat_ind++] = (data);
  }
  
}

void cprint(double c) {
  idouble cc;
  cc.d = c;

  ifloat real, imag;
  real.u32 = (cc.u64 >> 32);
  imag.u32 = (cc.u64 & 0xffffffff);
  std::cout << " (" << real.f << ")+(" << imag.f << "i) ";
}

// file -> gates vector
void parse(int num_qubits, int num_states, int max_num_gates, std::string filename, 
  std::vector<TYPE>& states, std::vector<Gate*>& gates, uint32_t& n_mat, uint32_t& n_ind) {
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
  for (int i = 1; i < num_states; i++)
     states[i] = 0;
  //   states[2 * i + 1] = 0;
  // }
  // states[0] = 1;
  states[0] = 0.0078125;

  std::string gate_st;
  std::string qubits_st;
  
  while (infile >> gate_st >> qubits_st) {
    std::vector<int> qubit_pos;
    int start = 0, end = 0;
    while(end < qubits_st.size()) {
      if(isdigit(qubits_st[end])) {
        end++;
      } else {
        if (start != end)
          qubit_pos.push_back(stoi(qubits_st.substr(start, end)));
        end++;
        start = end;
      }
    }
    Gate* g = generate_gates(gate_st, qubit_pos, gates);
    gates.push_back(g);
    n_ind += g->qubits.size();
    n_mat += g->matrix.size();
  }



  // std::cout << "\n\nAfter Parsing" << "\n";
  // for (uint32_t i = 0; i < num_qubits; i++) {
  //   std::cout << "Qubit-" << i << "  ";
  //   for (uint32_t g = 0; g < ng_per_qubit[i]; g++) {
  //     std::cout << gate_info[i][g] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // std::cout << "\n\n";
}

int main(int argc, char *argv[]) {
  // parse command arguments
  parse_args(argc, argv);

  // open device connection
  std::cout << "open device connection" << std::endl;
  RT_CHECK(vx_dev_open(&device));

  // std::srand(50);
  // Todo: rel add
  std::string filename = "/home/jm/vortex/tests/regression/qsim/temp.qasm";
  std::ifstream infile(filename);

  // first 4 lines
  std::string temp;
  for (int i = 0; i < 8; i++) { 
    infile >> temp;
  }
  uint32_t num_qubits = temp[temp.length() - 3] - '0'; 
  uint32_t max_num_gates = 30;
  uint32_t num_states = 1<<num_qubits;
  uint32_t num_indexes = 0;
  uint32_t num_matrices = 0;


  // generate source data
  // store state vector of each qubit
  std::vector<TYPE> h_A(num_states);
  // std::vector<TYPE> h_C(num_states);

  // Store gate array such that 
  // store gate array on each qubit
  std::vector<Gate*> GateArray;
  // // result state vector
  // std::vector<TYPE> h_C(2 * num_qubits);
  std::cout << temp << std::endl;
  std::cout << num_qubits << std::endl;
  
  parse(num_qubits, num_states, max_num_gates, filename, h_A, GateArray, num_matrices, num_indexes);

  uint32_t num_gates = GateArray.size();

  std::vector<int> h_num_q(num_gates + 1);
  std::vector<int> h_ind(num_indexes);
  std::vector<TYPE> h_mat(num_matrices);

  h_num_q[num_gates] = -1;

  convert_gates_to_data(GateArray, h_num_q, h_ind, h_mat);

  int j = 0;
  int k = 0;
  for(int i = 0; i < num_gates; ++i) {
    int q_num = h_num_q[i];
    std::cout << q_num << std::endl;
    
    for(int tmp = j;j < tmp + q_num; ++j)
      std::cout << h_ind[j] << ",";
    std::cout << std::endl;

    for(int tmp = k;k < tmp + (1<<q_num) * (1<<q_num);++k)
      // std::cout << h_mat[k] << "_";
      cprint(h_mat[k]);
    std::cout << std::endl;
  }

  uint32_t states_buf_size = num_states * sizeof(TYPE);
  uint32_t nqubit_buf_size = (num_gates + 1) * sizeof(int);
  uint32_t ind_buf_size = num_indexes * sizeof(int);
  uint32_t mat_buf_size = num_matrices * sizeof(TYPE);

  std::cout << "data type: " << Comparator<TYPE>::type_str() << std::endl;

  uint64_t num_cores, num_warps, num_threads;
  RT_CHECK(vx_dev_caps(device, VX_CAPS_NUM_CORES, &num_cores));
  RT_CHECK(vx_dev_caps(device, VX_CAPS_NUM_WARPS, &num_warps));
  RT_CHECK(vx_dev_caps(device, VX_CAPS_NUM_THREADS, &num_threads));

  std::cout << "num qubits " << num_qubits << std::endl;
  std::cout << "num gates " << num_gates << std::endl;
  std::cout << "num states " << num_states << std::endl;
  std::cout << "num cores " << num_cores << std::endl;
  std::cout << "num warps " << num_warps << std::endl;
  std::cout << "num threads " << num_threads << std::endl;

  kernel_arg.num_tasks = num_cores * num_warps * num_threads;
  kernel_arg.max_num_gates = max_num_gates;
  kernel_arg.num_qubits = num_qubits;
  kernel_arg.num_local_qubits = num_qubits - int(log(num_cores) / log(2));
  kernel_arg.num_states = num_states;

  std::cout << "allocate device memory" << std::endl;
  RT_CHECK(vx_dev_caps(device, VX_CAPS_LOCAL_MEM_ADDR, &kernel_arg.local_addr));
  RT_CHECK(vx_mem_alloc(device, states_buf_size, VX_MEM_READ_WRITE, &state_buffer));
  RT_CHECK(vx_mem_address(state_buffer, &kernel_arg.states_addr));
  RT_CHECK(vx_mem_alloc(device, nqubit_buf_size, VX_MEM_READ, &nqubit_buffer));
  RT_CHECK(vx_mem_address(nqubit_buffer, &kernel_arg.gate_num_qubits_addr));
  RT_CHECK(vx_mem_alloc(device, ind_buf_size, VX_MEM_READ_WRITE, &ind_buffer));
  RT_CHECK(vx_mem_address(ind_buffer, &kernel_arg.gate_qubits_addr));
  RT_CHECK(vx_mem_alloc(device, mat_buf_size, VX_MEM_READ, &matrix_buffer));
  RT_CHECK(vx_mem_address(matrix_buffer, &kernel_arg.gate_matrix_addr));
  // RT_CHECK(vx_mem_alloc(device, states_buf_size, VX_MEM_WRITE, &out_state_buffer));
  // RT_CHECK(vx_mem_address(out_state_buffer, &kernel_arg.out_states_addr));

  std::cout << "local_addr=0x" << std::hex << kernel_arg.local_addr << std::endl;
  std::cout << "state_addr=0x" << std::hex << kernel_arg.states_addr << std::endl;
  std::cout << "nqubits_addr=0x" << std::hex << kernel_arg.gate_num_qubits_addr << std::endl;
  std::cout << "ind_addr=0x" << std::hex << kernel_arg.gate_qubits_addr << std::endl;
  std::cout << "mat_addr=0x" << std::hex << kernel_arg.gate_matrix_addr << std::endl;
  // std::cout << "out_state_addr=0x" << std::hex << kernel_arg.out_states_addr << std::endl;

  // upload data
  {
    std::cout << "upload state buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(state_buffer, h_A.data(), 0, states_buf_size));
  }
  {
    std::cout << "upload q num buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(nqubit_buffer, h_num_q.data(), 0, nqubit_buf_size));
  }
  {
    std::cout << "upload ind buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(ind_buffer, h_ind.data(), 0, ind_buf_size));
  }
  {
    std::cout << "upload gate matrix buffer" << std::endl;
    RT_CHECK(vx_copy_to_dev(matrix_buffer, h_mat.data(), 0, mat_buf_size));
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

  // // download destination buffer
  std::cout << "download destination buffer" << std::endl;
  RT_CHECK(vx_copy_from_dev(h_A.data(), state_buffer, 0, states_buf_size));


  std::cout << "\n\n";
  for (uint32_t i = 0; i < num_states; i++) {
    std::cout << "ind " << i << " " ;
    cprint(h_A[i]);
    std::cout << std::endl;
  }
  std::cout << "\n\n";


  // cleanup
  std::cout << "cleanup" << std::endl;
  cleanup();

  std::cout << "PASSED!" << std::endl;

  return 0;
}
