OPENQASM 2.0;
include "temp.inc";
qreg q[5];
creg c[5];

h q[0];
cx q[0],q[4];