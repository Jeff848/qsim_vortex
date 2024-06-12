#include <stdint.h>
#include <vx_intrinsics.h>
#include <vx_spawn.h>
#include <vx_print.h>
#include "common.h"

inline char is_log2(uint32_t x) {
    return ((x & (x-1)) == 0);
}

inline complex mul(complex c1, complex c2) {
  complex temp;
  temp.real = c1.real * c2.real - c1.imag * c2.imag;
  temp.imag = c1.real * c2.imag - c1.imag * c2.real;
  return (temp);
}

inline complex add(complex c1, complex c2) {
  complex temp;
  temp.real = c1.real + c2.real;
  temp.imag = c1.imag + c2.imag;
  return (temp);
}

// task_id: qubit number (qidx)
void kernel_body(uint32_t task_id, kernel_arg_t* __UNIFORM__ arg) {
	auto A = reinterpret_cast<TYPE*>(arg->A_addr);
	auto B = reinterpret_cast<TYPE*>(arg->B_addr);
	auto C = reinterpret_cast<TYPE*>(arg->C_addr);
  auto nq = arg->nq;



  auto size = arg->nst;
  for (int i = 0; i < arg->max_num_gates; i++) {
    for (int j = 0; j < size; j++) {
      TYPE rowsum = {0, 0};
      for (int k = 0; k < size; k++) {
        TYPE temp = mul(B[i * size * size + j * size + k], A[k]);
        rowsum = add(temp, rowsum);
      }
      C[j] = rowsum;
    }

    // for (int i = 0; i < size; i++) {
    //   TYPE k = C[i];
    //   A[i] = k;
    // }
  }
}

int main() {
	kernel_arg_t* arg = (kernel_arg_t*)csr_read(VX_CSR_MSCRATCH);
	vx_spawn_tasks(arg->num_tasks, (vx_spawn_tasks_cb)kernel_body, arg);
	return 0;
}