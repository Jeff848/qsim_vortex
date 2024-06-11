#include <stdint.h>
#include <vx_intrinsics.h>
#include <vx_spawn.h>
#include <vx_print.h>
#include "common.h"

inline char is_log2(uint32_t x) {
    return ((x & (x-1)) == 0);
}

// task_id: qubit number (qidx)
//For actual implementation should allocate at least 2 qubits per node
//Because we want to work with at least 2 qubit gates
//Gates array should have matrix size included
void kernel_swap(int group1, int group2) {
    auto num_cores = vx_num_cores();
	auto num_warps = vx_num_warps();
	auto num_threads = vx_num_threads();

	auto cid = vx_core_id();
	auto wid = vx_warp_id();
	auto tid = vx_thread_id();

    // memory fence
	vx_fence();

	// global barrier
	vx_barrier(0x80000000, num_cores);
}

//bit twiddling
uint32_t get_bit(int n, int t) {
    return (n >> t) & 1;
}

uint32_t set_bits(int n, int* qubits, int num_qubits, int v) {
    for (int q = 0; q < num_qubits; ++q) {
        int b = get_bit(v, q) << qubits[q];
        n = (n & (~b)) | b;
    }
    return n;
}

uint32_t insert_bits(int n, int* qubits, int num_qubits, int b) {
    for (int q = 0; q < num_qubits; ++q) {
        int t = qubits[q];
        int l = (n>>t)<<(t+1);
        int m = b << t;
        int r = n & ((1<<t) - 1);
        n = l | m | r;
    }
    return n;
}

//separate kernels for doing matmul and doing gate fusion
void kernel_body(uint32_t task_id, kernel_arg_t* __UNIFORM__ arg) {
    auto local_ptr = reinterpret_cast<TYPE*>(arg->local_addr);
    auto states = reinterpret_cast<TYPE*>(arg->states_addr);
	auto num_qubits_arr = reinterpret_cast<int*>(arg->gate_num_qubits_addr);
    auto gate_indexes_arr = reinterpret_cast<int*>(arg->gate_qubits_addr);
    auto gate_matrix_arr = reinterpret_cast<TYPE*>(arg->gate_matrix_addr);

    auto num_cores = vx_num_cores();
	auto num_warps = vx_num_warps();
	auto num_threads = vx_num_threads();

	auto cid = vx_core_id();
	auto wid = vx_warp_id();
	auto tid = vx_thread_id();
	// auto C = reinterpret_cast<TYPE*>(arg->C_addr);
    
    auto nq = arg->num_qubits;
    auto nstates = arg->num_states;

    // uint32_t st_b = 4 * arg->max_num_gates * task_id;       // starting index of gate array
    // uint32_t st_a = 2 * task_id;                            // starting index of state vector

    // initial vector of qubit qidx

    // local indexes
    
    // global indexes

    // swap

    // g: gate number
    uint32_t q = 0;
    uint32_t m = 0;

    uint32_t global_ind = cid;
    uint32_t local_ind = wid * num_threads + tid;

    vx_printf("global=%d\n", global_ind);
    vx_printf("local=%d\n", local_ind);


    uint32_t states_per_core = nstates / num_cores;
    auto local_states = local_ptr;
    vx_printf("local_states=%x\n", local_states);

    //Per core mapping of original order qubit to swap order


    for (uint32_t g = 0; g < (arg->max_num_gates); g++) {
        // end of gate array
        if (num_qubits_arr[g] == -1) {
            break;
        }
        uint32_t num_qubits = num_qubits_arr[g];
        uint32_t num_data = (1 << num_qubits);
        
        //swap if global index
        //global barrier and swap
        
        // matrix - vector multiplication
        

        //Get number of gate operations on core
        uint32_t num_local_applications = states_per_core / num_data;
        uint32_t num_local_per_thread = num_local_applications / (num_warps * num_threads);
        vx_printf("num_local_applications %d\n", num_local_applications);
        //if num_local_per_thread is less than 0
        if (num_local_per_thread > 0) {
            for(uint32_t i = local_ind * num_local_per_thread; 
                i < local_ind * num_local_per_thread + num_local_per_thread; ++i) {

                //sync across 
                //map i to permutation
                
                int z = insert_bits(i, gate_indexes_arr + q, num_qubits, 0);

                for(uint32_t j = 0; j < num_data; ++j) {
                    int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                    local_states[i * num_data + j] = states[global_ind * states_per_core + ind];
                }

                vx_fence();
                // local barrier
                vx_barrier(0, num_warps);

                
                for(uint32_t j = 0; j < num_data; ++j) {
                    int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                    TYPE temp = 0;
                    for(uint32_t l = 0; l < num_data; ++l) {
                        temp += gate_matrix_arr[m + j * num_data + l] * local_states[i * num_data + l];
                    }
                    states[global_ind * states_per_core + ind] = temp;
                }
                
                vx_fence();
                // local barrier
                vx_barrier(0, num_warps);
            }
        } else {
            
            //Assume that we assign each warp to a gate application 
            uint32_t num_thread_per_app = (num_threads * num_warps) / num_local_applications;
            uint32_t num_data_per_thread = num_data / num_thread_per_app;
            uint32_t local_application_ind = local_ind % num_local_applications; 
            uint32_t app_ind = local_ind / num_local_applications;

            int z = insert_bits(local_application_ind, gate_indexes_arr + q, num_qubits, 0);
            vx_printf("z: %d\n", z);
            vx_printf("local_application_ind %d\n", local_application_ind);
            vx_printf("app_ind %d\n", app_ind);
            for(uint32_t j = num_data_per_thread * app_ind; j < num_data_per_thread * app_ind + num_data_per_thread; ++j) {
                int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                local_states[local_application_ind * num_data + j] = states[global_ind * states_per_core + ind];
            }
            
            vx_fence();
            // local barrier
            vx_barrier(0, num_warps);
            
            for(uint32_t j = num_data_per_thread * app_ind; j < num_data_per_thread * app_ind + num_data_per_thread; ++j) {
                int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                TYPE temp = 0;
                for(uint32_t l = 0; l < num_data; ++l) {
                    temp += gate_matrix_arr[m + j * num_data + l] * local_states[local_application_ind * num_data + l];
                }
                states[global_ind * states_per_core + ind] = temp;
            }

            vx_fence();
            // local barrier
            vx_barrier(0, num_warps); 
        }

        q += num_qubits; 
        m += num_data * num_data;

        //sync across warps
        // memory fence
        // vx_fence();
        // // local barrier
        // vx_barrier(0, num_warps);
    }

    vx_printf("task=%d\n", task_id);
}

int main() {
	kernel_arg_t* arg = (kernel_arg_t*)csr_read(VX_CSR_MSCRATCH);
    
    
	vx_spawn_tasks(arg->num_tasks, (vx_spawn_tasks_cb)kernel_body, arg);
	return 0;
}
