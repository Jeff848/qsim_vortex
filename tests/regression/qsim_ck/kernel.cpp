#include <stdint.h>
#include <vx_intrinsics.h>
#include <vx_spawn.h>
#include <vx_print.h>
#include "common.h"

inline char is_log2(uint32_t x) {
    return ((x & (x-1)) == 0);
}

// task_id: qubit number (qidx)
void kernel_body(uint32_t task_id, kernel_arg_t* __UNIFORM__ arg) {
	auto A = reinterpret_cast<TYPE*>(arg->A_addr);
	auto B = reinterpret_cast<TYPE*>(arg->B_addr);
	auto C = reinterpret_cast<TYPE*>(arg->C_addr);
    auto nq = arg->nq;

    uint32_t st_b = 4 * arg->max_num_gates * task_id;       // starting index of gate array
    uint32_t st_a = 2 * task_id;                            // starting index of state vector

    // initial vector of qubit qidx
    TYPE v1 = A[st_a];
    TYPE v2 = A[st_a + 1];

    // g: gate number
    for (uint32_t g = 0; g < (arg->max_num_gates); g++) {
        // end of gate array
        if (B[st_b + 4 * g] == 100) {
            break;
        }

        // matrix - vector multiplication
        // applying gates
        TYPE nv1 = B[st_b + 4*g] * v1 + B[st_b + 4*g + 1] * v2;
        TYPE nv2 = B[st_b + 4*g + 2] * v1 + B[st_b + 4*g + 3] * v2;

        v1 = nv1;
        v2 = nv2;
    }

    C[st_a] = v1;
    C[st_a + 1] = v2;
    // vx_printf("task=%d\n", task_id);
}

int main() {
	kernel_arg_t* arg = (kernel_arg_t*)csr_read(VX_CSR_MSCRATCH);
	vx_spawn_tasks(arg->num_tasks, (vx_spawn_tasks_cb)kernel_body, arg);
	return 0;
}
