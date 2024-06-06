#ifndef _COMMON_H_
#define _COMMON_H_

#ifndef TYPE
#define TYPE float
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



#endif
