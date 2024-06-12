#ifndef _COMMON_H_
#define _COMMON_H_

// #ifndef TYPE
// #define TYPE float
// #endif

typedef struct {
  uint32_t num_tasks;
  uint32_t max_num_gates;
  uint32_t nst;
  uint32_t nq;
  uint64_t A_addr;
  uint64_t B_addr;
  uint64_t C_addr;  
} kernel_arg_t;

typedef struct {
  double real;
  double imag;
} complex;

#ifndef TYPE
#define TYPE complex
#endif

#endif