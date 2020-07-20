`include "VX_define.vh"

module VX_writeback #(
    parameter CORE_ID = 0
) (
    input wire      clk,
    input wire      reset,

    // inputs
    VX_commit_if    alu_commit_if,
    VX_commit_if    branch_commit_if,
    VX_commit_if    lsu_commit_if,  
    VX_commit_if    mul_commit_if,    
    VX_commit_if    csr_commit_if,

    // outputs
    VX_wb_if        writeback_if
);

    wire br_valid  = (| branch_commit_if.valid) && (branch_commit_if.wb != `WB_NO);
    wire lsu_valid = (| lsu_commit_if.valid) && (lsu_commit_if.wb != `WB_NO);
    wire mul_valid = (| mul_commit_if.valid) && (mul_commit_if.wb != `WB_NO);
    wire alu_valid = (| alu_commit_if.valid) && (alu_commit_if.wb != `WB_NO);
    wire csr_valid = (| csr_commit_if.valid) && (csr_commit_if.wb != `WB_NO);

    VX_wb_if writeback_tmp_if();    

    assign writeback_tmp_if.valid =  br_valid ? branch_commit_if.valid :
                                    lsu_valid ? lsu_commit_if.valid :
                                    mul_valid ? mul_commit_if.valid :             
                                    alu_valid ? alu_commit_if.valid :                            
                                    csr_valid ? csr_commit_if.valid :                                                 
                                                0;     

    assign writeback_tmp_if.warp_num = br_valid ? branch_commit_if.warp_num :
                                    lsu_valid ? lsu_commit_if.warp_num :
                                    mul_valid ? mul_commit_if.warp_num :   
                                    alu_valid ? alu_commit_if.warp_num :                            
                                    csr_valid ? csr_commit_if.warp_num :                           
                                    
                                                0; 

    assign writeback_tmp_if.data =   br_valid ? branch_commit_if.data :
                                    lsu_valid ? lsu_commit_if.data :
                                    mul_valid ? mul_commit_if.data :                           
                                    alu_valid ? alu_commit_if.data :                            
                                    csr_valid ? csr_commit_if.data :                                                               
                                                0;

    assign writeback_tmp_if.rd =     br_valid ? branch_commit_if.rd :
                                    lsu_valid ? lsu_commit_if.rd :
                                    mul_valid ? mul_commit_if.rd :                           
                                    alu_valid ? alu_commit_if.rd :                            
                                    csr_valid ? csr_commit_if.rd :                                                               
                                                0;

    wire stall = ~writeback_if.ready && (| writeback_if.valid);

    VX_generic_register #(
        .N(`NUM_THREADS + `NW_BITS + `NR_BITS + (`NUM_THREADS * 32))
    ) wb_reg (
        .clk   (clk),
        .reset (reset),
        .stall (stall),
        .flush (0),
        .in    ({writeback_tmp_if.valid, writeback_tmp_if.warp_num, writeback_tmp_if.rd, writeback_tmp_if.data}),
        .out   ({writeback_if.valid,     writeback_if.warp_num,     writeback_if.rd,     writeback_if.data})
    );

    assign branch_commit_if.ready = !stall;    
    assign lsu_commit_if.ready    = !stall && !br_valid;    
    assign mul_commit_if.ready    = !stall && !br_valid && !lsu_valid;
    assign alu_commit_if.ready    = !stall && !br_valid && !lsu_valid && !mul_valid;    
    assign csr_commit_if.ready    = !stall && !br_valid && !lsu_valid && !mul_valid && !alu_valid;    
    
    // special workaround to control RISC-V benchmarks termination on Verilator
    reg [31:0] last_data_wb /* verilator public */;
    always @(posedge clk) begin
        if ((| writeback_tmp_if.valid) && ~stall && (writeback_tmp_if.rd == 28)) begin
            last_data_wb <= writeback_tmp_if.data[0];
        end
    end

endmodule






