#ifndef _COMMON_H_
#define _COMMON_H_

#ifndef TYPE
#define TYPE double
#endif

typedef struct {
  uint32_t num_tasks;
  uint32_t max_num_gates;
  uint32_t num_qubits;
  uint32_t num_states;
  uint64_t local_addr;
  uint64_t states_addr;
  uint64_t gate_num_qubits_addr;
  uint64_t gate_qubits_addr; 
  uint64_t gate_matrix_addr;
} kernel_arg_t;

typedef union {
  double d;
  int64_t i64;
  uint64_t u64;
} idouble;

typedef union {
  float f;
  int32_t i32;
  uint32_t u32;
} ifloat;

#endif
