#include <stdint.h>
#include <vx_intrinsics.h>
#include <vx_spawn.h>
#include <vx_print.h>
#include "common.h"

inline char is_log2(uint32_t x) {
    return ((x & (x-1)) == 0);
}



//bit twiddling
uint32_t flip_bit(int n, int t) {
    return n ^ (1<<t);
}

uint32_t flip_bits(int n, int* qubits, int num_qubits) {
    for (int q = 0; q < num_qubits; ++q) {
        n = flip_bit(n, qubits[q]);
    }
    return n;
}


uint32_t get_bitmask(int* qubits, int num_qubits) {
    return flip_bits(0, qubits, num_qubits);
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

uint32_t insert_bit(int n, int t, int b) {
    int l = (n>>t)<<(t+1);
    int m = b << t;
    int r = n & ((1<<t) - 1);
    n = l | m | r;
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

uint32_t swap_bits(int n, int p1, int p2) {
   //left-shift 1 p1 and p2 times 
   //and using XOR
  if (((n & (1 << p1)) >> p1) ^ ((n & (1 << p2)) >> p2)) {
    n ^= 1 << p1;
    n ^= 1 << p2;
  }
  return n;
}

uint32_t swap_bits_arr(int n, int* qubits, int* swap_qubits, int num_qubits) {
    for (int q = 0; q < num_qubits; ++q) {
        if (qubits[q] != swap_qubits[q]) {
            n = swap_bits(n, qubits[q], swap_qubits[q]);
        }
    }
    return n;
}



uint32_t set_swap_bits(int n, int* qubits, int* swap_qubits, int num_qubits, int v) {
    int i = 0;
    for (int q = 0; q < num_qubits; ++q) {
        if (qubits[q] != swap_qubits[q]) {
            int swap;
            if (qubits[q] < swap_qubits[q]) {
                swap = qubits[q];
            } else {
                swap = swap_qubits[q];
            }
            int b = get_bit(v, i) << swap;
            n = (n & (~b)) | b;
            i = i+1;
        }
    }
    return n;
}

void swap_op(TYPE* states, int* swap_qubits, int* qubit_indexes, int num_local_ind, TYPE* local_addr, int global_ind, int target_core, int num_qubits) {
    //get bit vector set with the swap qubits that aren't equal to the qubit indexes
    int j = set_swap_bits(0, qubit_indexes, swap_qubits, num_qubits, target_core); 
    vx_printf("swap bit vector %d\n", j);
    for(int i = (target_core << num_local_ind); i < (target_core << num_local_ind) + (1<<num_local_ind); ++i) {

        if((j & i) ^ j) { //value should be swapped
            int r = swap_bits_arr(i, qubit_indexes, swap_qubits, num_qubits);
            vx_printf("i r %d, %d\n", i, r);
            // do swap
            TYPE temp;
            temp = states[i];
            states[i] = states[r];
            states[r] = temp;
            // local_addr[global_ind * states_per_core + r] = states[target_core * states_per_core + i];
        } 
    }    
}

// task_id: qubit number (qidx)
//For actual implementation should allocate at least 2 qubits per node
//Because we want to work with at least 2 qubit gates
//Gates array should have matrix size included
void do_swap(int* qubit_indexes, int num_qubits, int num_local_ind, int* swap_targets, TYPE* states, TYPE* local_addr, int global_ind, int num_states) {
    auto num_cores = vx_num_cores();
	auto num_warps = vx_num_warps();
	auto num_threads = vx_num_threads();

	auto cid = vx_core_id();
	auto wid = vx_warp_id();
	auto tid = vx_thread_id();

    for(int i = 1; i < num_cores; ++i) {
        vx_printf("sv pair: %d, %d\n", i, i ^ global_ind);
        int pair_core = i ^ global_ind;
        swap_op(states, swap_targets, qubit_indexes, num_local_ind, local_addr, global_ind, pair_core, num_qubits);

        // memory fence
        vx_fence();

        // global barrier
        vx_barrier(0x80000000, num_cores);

        
    }

    for(int q = 0; q < num_qubits; ++q) {
        int temp;
        temp = qubit_indexes[q];
        qubit_indexes[q] = swap_targets[q];
        swap_targets[q] = temp;
    }

}

void bit_swap(int* qubit_indexes, int num_qubits, int num_local_ind, int* local_addr) {
    auto num_cores = vx_num_cores();
	auto num_warps = vx_num_warps();
	auto num_threads = vx_num_threads();

	auto cid = vx_core_id();
	auto wid = vx_warp_id();
	auto tid = vx_thread_id();


    //Calculate earliest non-target qubit
    uint32_t b = get_bitmask(qubit_indexes, num_qubits);
    uint32_t q = 0;
    while (get_bit(b, q) == 1)
        q += 1;

    //Which targets need switching
    int* swap_targets = local_addr;
    for(int i = 0; i < num_qubits; ++i) {
        if (qubit_indexes[i] < num_local_ind)
            swap_targets[i] = qubit_indexes[i];
        else {

            swap_targets[i] = q;
            q += 1;
            while (get_bit(b, q) == 1)
                q += 1;
        }
    }
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

    auto local_swap_arr = reinterpret_cast<int*>(arg->local_addr + nstates + 1000);

    // g: gate number
    uint32_t q = 0;
    uint32_t m = 0;

    uint32_t global_ind = cid;
    uint32_t local_ind = wid * num_threads + tid;
    uint32_t nlocalq = arg -> num_local_qubits;

    vx_printf("global=%d\n", global_ind);
    vx_printf("local=%d\n", local_ind);
    vx_printf("localq=%d\n", nlocalq);

    uint32_t states_per_core = nstates / num_cores;
    auto local_states = local_ptr;
    vx_printf("local_states=%x\n", local_states);


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


        bit_swap(gate_indexes_arr + q, num_qubits, nlocalq, local_swap_arr);

        for(int i = 0; i < num_qubits; ++i) {
            vx_printf("swap %d ", local_swap_arr[i]);
        }
        vx_printf("\n");

        for(int i = 0; i < num_qubits; ++i) {
            vx_printf("qub : %d, ", (gate_indexes_arr + q)[i]);
        }
        vx_printf("\n");

        do_swap(gate_indexes_arr + q, num_qubits, nlocalq, local_swap_arr, states, local_states, global_ind, nstates);

        for(int i = 0; i < num_qubits; ++i) {
            vx_printf("qub : %d, ", (gate_indexes_arr + q)[i]);
        }
        vx_printf("\n");
        

        if (num_local_per_thread > 0) {
            for(uint32_t i = local_ind * num_local_per_thread; 
                i < local_ind * num_local_per_thread + num_local_per_thread; ++i) {

                
                int z = insert_bits(i, gate_indexes_arr + q, num_qubits, 0);
  
                for(uint32_t j = 0; j < num_data; ++j) {
                    int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                    local_states[i * num_data + j] = states[global_ind * states_per_core + ind];
                }

                
                for(uint32_t j = 0; j < num_data; ++j) {
                    int ind = set_bits(z, gate_indexes_arr + q, num_qubits, j);
                    TYPE temp = 0;
                    for(uint32_t l = 0; l < num_data; ++l) {
                        TYPE temp_mul = vx_cmul(gate_matrix_arr[m + j * num_data + l], local_states[i * num_data + l]);
                        temp = vx_cadd(temp_mul, temp);
                    }
                    states[global_ind * states_per_core + ind] = temp;
                }
                
               
            }
        } else {
            
            //Assume that we assign each warp to a gate application 
            uint32_t num_thread_per_app = (num_threads * num_warps) / num_local_applications;
            uint32_t num_data_per_thread = num_data / num_thread_per_app;
            uint32_t local_application_ind = local_ind / (states_per_core / num_local_applications); 
            uint32_t app_ind = local_ind % (states_per_core / num_local_applications);

            int z = insert_bits(local_application_ind, gate_indexes_arr + q, num_qubits, 0);
            // vx_printf("z: %d\n", z);
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
                    // temp += gate_matrix_arr[m + j * num_data + l] * local_states[local_application_ind * num_data + l];
                    TYPE temp_mul = vx_cmul(gate_matrix_arr[m + j * num_data + l], local_states[local_application_ind * num_data + l]);
                    temp = vx_cadd(temp_mul, temp);
                }
                states[global_ind * states_per_core + ind] = temp;
            }
        }

        do_swap(gate_indexes_arr + q, num_qubits, nlocalq, local_swap_arr, states, local_states, global_ind, nstates);
        for(int i = 0; i < num_qubits; ++i) {
            vx_printf("swap %d ", local_swap_arr[i]);
        }
        vx_printf("\n");

        for(int i = 0; i < num_qubits; ++i) {
            vx_printf("qub : %d, ", (gate_indexes_arr + q)[i]);
        }

        q += num_qubits; 
        m += num_data * num_data;

        // do_swap();
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
