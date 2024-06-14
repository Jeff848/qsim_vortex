#ifndef GATES
#define GATES

#include <vector>
#include "common.h"

struct gate_data {
  int num_qubits;
  int* qubits;
  TYPE* matrix;
}; 

//consists of qubit indexes and a generate_matrix function
//associated matrix
class Gate {
public:
    std::vector<int> qubits;
    std::vector<TYPE> matrix; //store data in column-order
    // std::string gate_name;
    virtual void initialize(std::vector<int> qubit_ordering) = 0;
    virtual void initialize_matrix() = 0;



    Gate(int mat_size, int def_value, int num_qubits) : 
        qubits(num_qubits),
        matrix(mat_size, def_value)
    {}
    
    // void merge_gates(Gate& a, Gate& b) {
    //     //qubits is union of a and b
    //     //matrix is fusion of a and b
    //     //If fusion is large enough then maybe send off to the GPU?
    // }
    
};

class OneQubitGate : public Gate { 
public:
    OneQubitGate() : Gate(4, 0, 1) {
    }
    void initialize(std::vector<int> qubit_ordering) {
        qubits = qubit_ordering;
        //Verify that qubits is len 1
        this->initialize_matrix();
    }

    virtual void initialize_matrix() = 0;
};

class TwoQubitGate : public Gate {
public:
    TwoQubitGate() : Gate(16, 0, 2) {
    }

    void initialize(std::vector<int> qubit_ordering) {
        qubits = qubit_ordering;
        
        //Verify that qubits is len 2
        if(qubits[0] > qubits[1]) {//Sort the qubit order
            int temp = qubits[0];
            qubits[0] = qubits[1];
            qubits[1] = temp;
            this->initialize_matrix_reverse();
        } else {
            this->initialize_matrix();
        }
        
    }

    virtual void initialize_matrix() = 0;
    virtual void initialize_matrix_reverse() = 0;
};

class HGate: public OneQubitGate {
public:
    HGate(): OneQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "H";
        ifloat f1, f2, f3;
        f1.f = 0.70710678118;
        f2.f = -0.70710678118;
        f3.f = 0;

        idouble d1, d2;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f3.u32;
        d2.u64 = ((uint64_t)(f2.u32) << 32) | f3.u32;

        this->matrix[0] = d1.d;
        this->matrix[1] = d1.d;
        this->matrix[2] = d1.d;
        this->matrix[3] = d2.d;
    }
};

class XGate: public OneQubitGate {
public:
    XGate(): OneQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "X";
        ifloat f1, f2;
        f1.f = 1;
        f2.f = 0;

        idouble d1;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f2.u32;

        this->matrix[1] = d1.d;
        this->matrix[2] = d1.d;
    }
};

class ZGate: public OneQubitGate {
public:
    ZGate(): OneQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "Z";
        ifloat f1, f2, f3;
        f1.f = 1;
        f2.f = -1;
        f3.f = 0;

        idouble d1, d2;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f3.u32;
        d2.u64 = ((uint64_t)(f2.u32) << 32) | f3.u32;

        this->matrix[1] = d1.d;
        this->matrix[2] = d2.d;
    }
};


class TGate: public OneQubitGate {
public:
    TGate(): OneQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "T";
        ifloat f1, f2, f3;
        f1.f = 1;
        f2.f = 0.70710678118;
        f3.f = 0;

        idouble d1, d2, d3;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f3.u32;
        d2.u64 = ((uint64_t)(f2.u32) << 32) | f2.u32;

        this->matrix[0] = d1.d;
        this->matrix[3] = d2.d;
    }
};

class TdgGate: public OneQubitGate {
public:
    TdgGate(): OneQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "Tdg";
        ifloat f1, f2, f3, f4;
        f1.f = 1;
        f2.f = 0.70710678118;
        f3.f = 0;
        f4.f = -0.70710678118;

        idouble d1, d2, d3;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f3.u32;
        d2.u64 = ((uint64_t)(f2.u32) << 32) | f4.u32;

        this->matrix[0] = d1.d;
        this->matrix[3] = d2.d;
    }
};

class CxGate: public TwoQubitGate {
public:
    CxGate(): TwoQubitGate()
    {}

    void initialize_matrix() {
        // this->gate_name = "CX";
        ifloat f1, f2;
        f1.f = 1;
        f2.f = 0;

        idouble d1;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f2.u32;

        this->matrix[0] = d1.d;
        this->matrix[6] = d1.d;
        this->matrix[11] = d1.d;
        this->matrix[13] = d1.d;
    }
    void initialize_matrix_reverse() {
        // this->gate_name = "CX";

        ifloat f1, f2;
        f1.f = 1;
        f2.f = 0;

        idouble d1;
        d1.u64 = ((uint64_t)(f1.u32) << 32) | f2.u32;

        this->matrix[0] = d1.d;
        this->matrix[5] = d1.d;
        this->matrix[11] = d1.d;
        this->matrix[14] = d1.d;
    }
};

// class GateWrapper {
// public:
//     HGate h;
//     XGate x;
//     TGate t;
//     TdgGate tdg;
//     CxGate cx;
 
//     void generate_gates(std::vector<TYPE>& qubit_ordering) {
        
//     }
// }

#endif

// #ifndef GATES
// #define GATES

// #include <vector>
// #include "common.h"

// struct gate_data {
//   int num_qubits;
//   int* qubits;
//   TYPE* matrix;
// }; 

// //consists of qubit indexes and a generate_matrix function
// //associated matrix
// class Gate {
// public:
//     std::vector<int> qubits;
//     std::vector<TYPE> matrix; //store data in column-order
//     // std::string gate_name;
//     virtual void initialize(std::vector<int> qubit_ordering) = 0;
//     virtual void initialize_matrix() = 0;



//     Gate(int mat_size, int def_value, int num_qubits) : 
//         qubits(num_qubits),
//         matrix(mat_size, def_value)
//     {}
    
//     // void merge_gates(Gate& a, Gate& b) {
//     //     //qubits is union of a and b
//     //     //matrix is fusion of a and b
//     //     //If fusion is large enough then maybe send off to the GPU?
//     // }
    
// };

// class OneQubitGate : public Gate { 
// public:
//     OneQubitGate() : Gate(4, 0, 1) {
//     }
//     void initialize(std::vector<int> qubit_ordering) {
//         qubits = qubit_ordering;
//         //Verify that qubits is len 1
//         this->initialize_matrix();
//     }

//     virtual void initialize_matrix() = 0;
// };

// class TwoQubitGate : public Gate {
// public:
//     TwoQubitGate() : Gate(16, 0, 2) {
//     }

//     void initialize(std::vector<int> qubit_ordering) {
//         qubits = qubit_ordering;
        
//         //Verify that qubits is len 2
//         if(qubits[0] > qubits[1]) {//Sort the qubit order
//             int temp = qubits[0];
//             qubits[0] = qubits[1];
//             qubits[1] = temp;
//             this->initialize_matrix_reverse();
//         } else {
//             this->initialize_matrix();
//         }
        
//     }

//     virtual void initialize_matrix() = 0;
//     virtual void initialize_matrix_reverse() = 0;
// };

// class HGate: public OneQubitGate {
// public:
//     HGate(): OneQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "H";
//         this->matrix[0] = 0.70710678118;
//         this->matrix[1] = 0.70710678118;
//         this->matrix[2] = 0.70710678118;
//         this->matrix[3] = -0.70710678118;
//     }
// };

// class XGate: public OneQubitGate {
// public:
//     XGate(): OneQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "X";
//         this->matrix[1] = 1;
//         this->matrix[2] = 1;
//     }
// };

// class ZGate: public OneQubitGate {
// public:
//     ZGate(): OneQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "Z";
//         this->matrix[0] = 1;
//         this->matrix[3] = -1;
//     }
// };


// class TGate: public OneQubitGate {
// public:
//     TGate(): OneQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "T";
//         this->matrix[0] = 1;
//         this->matrix[3] = -0.70710678118;
//     }
// };

// class TdgGate: public OneQubitGate {
// public:
//     TdgGate(): OneQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "Tdg";
//         this->matrix[0] = 0.70710678118;
//         this->matrix[1] = 0.70710678118;
//         this->matrix[2] = 0.70710678118;
//         this->matrix[3] = -0.70710678118;
//     }
// };

// class CxGate: public TwoQubitGate {
// public:
//     CxGate(): TwoQubitGate()
//     {}

//     void initialize_matrix() {
//         // this->gate_name = "CX";
//         this->matrix[0] = 1;
//         this->matrix[6] = 1;
//         this->matrix[11] = 1;
//         this->matrix[13] = 1;
//     }
//     void initialize_matrix_reverse() {
//         // this->gate_name = "CX";
//         this->matrix[0] = 1;
//         this->matrix[5] = 1;
//         this->matrix[11] = 1;
//         this->matrix[14] = 1;
//     }
// };

// // class GateWrapper {
// // public:
// //     HGate h;
// //     XGate x;
// //     TGate t;
// //     TdgGate tdg;
// //     CxGate cx;
 
// //     void generate_gates(std::vector<TYPE>& qubit_ordering) {
        
// //     }
// // }

// #endif