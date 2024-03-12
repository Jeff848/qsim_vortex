// Copyright © 2019-2023
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <future>
#include <list>
#include <chrono>

#include <vortex.h>
#include <malloc.h>
#include <utils.h>
#include <VX_config.h>
#include <VX_types.h>

#include <mem.h>
#include <util.h>
#include <processor.h>

#define RAM_PAGE_SIZE 4096

using namespace vortex;

///////////////////////////////////////////////////////////////////////////////

class vx_device {    
public:
    vx_device() 
        : ram_(0, RAM_PAGE_SIZE)
        , global_mem_(
            ALLOC_BASE_ADDR,
            ALLOC_MAX_ADDR - ALLOC_BASE_ADDR,
            RAM_PAGE_SIZE,
            CACHE_BLOCK_SIZE)
        , local_mem_(
            LMEM_BASE_ADDR,
            (1ull << LMEM_LOG_SIZE),
            RAM_PAGE_SIZE,
            1) 
    {
        processor_.attach_ram(&ram_);
    }

    ~vx_device() {    
        if (future_.valid()) {
            future_.wait();
        }
    }

    int mem_alloc(uint64_t size, int type, uint64_t* dev_addr) {
        if (type == VX_MEM_TYPE_GLOBAL) {
            return global_mem_.allocate(size, dev_addr);
        } else if (type == VX_MEM_TYPE_LOCAL) {
            return local_mem_.allocate(size, dev_addr);
        }
        return -1;
    }

    int mem_free(uint64_t dev_addr) {
        if (dev_addr >= LMEM_BASE_ADDR) {
            return local_mem_.release(dev_addr);
        } else {
            return global_mem_.release(dev_addr);
        }
    }

    int mem_info(int type, uint64_t* mem_free, uint64_t* mem_used) const {
        if (type == VX_MEM_TYPE_GLOBAL) {
            if (mem_free)
                *mem_free = global_mem_.free();
            if (mem_used)
                *mem_used = global_mem_.allocated();
        } else if (type == VX_MEM_TYPE_LOCAL) {
            if (mem_free)
                *mem_free = local_mem_.free();
            if (mem_used)
                *mem_free = local_mem_.allocated();
        } else {
            return -1;
        }
        return 0;
    }

    int upload(uint64_t dest_addr, const void* src, uint64_t size) {
        uint64_t asize = aligned_size(size, CACHE_BLOCK_SIZE);
        if (dest_addr + asize > GLOBAL_MEM_SIZE)
            return -1;

        /*printf("VXDRV: upload %ld bytes from 0x%lx:", size, uintptr_t((uint8_t*)src));
        for (int i = 0;  i < (asize / CACHE_BLOCK_SIZE); ++i) {
            printf("\n0x%08lx=", dest_addr + i * CACHE_BLOCK_SIZE);
            for (int j = 0;  j < CACHE_BLOCK_SIZE; ++j) {
                printf("%02x", *((uint8_t*)src + i * CACHE_BLOCK_SIZE + CACHE_BLOCK_SIZE - 1 - j));
            }
        }
        printf("\n");*/
        
        ram_.write((const uint8_t*)src, dest_addr, size);
        return 0;
    }

    int download(void* dest, uint64_t src_addr, uint64_t size) {
        uint64_t asize = aligned_size(size, CACHE_BLOCK_SIZE);
        if (src_addr + asize > GLOBAL_MEM_SIZE)
            return -1;

        ram_.read((uint8_t*)dest, src_addr, size);
        
        /*printf("VXDRV: download %ld bytes to 0x%lx:", size, uintptr_t((uint8_t*)dest));
        for (int i = 0;  i < (asize / CACHE_BLOCK_SIZE); ++i) {
            printf("\n0x%08lx=", src_addr + i * CACHE_BLOCK_SIZE);
            for (int j = 0;  j < CACHE_BLOCK_SIZE; ++j) {
                printf("%02x", *((uint8_t*)dest + i * CACHE_BLOCK_SIZE + CACHE_BLOCK_SIZE - 1 - j));
            }
        }
        printf("\n");*/
        
        return 0;
    }

    int start() {   
        // ensure prior run completed
        if (future_.valid()) {
            future_.wait();
        }
        // start new run
        future_ = std::async(std::launch::async, [&]{
            processor_.run();
        });
        return 0;
    }

    int wait(uint64_t timeout) {
        if (!future_.valid())
            return 0;
        uint64_t timeout_sec = timeout / 1000;
        std::chrono::seconds wait_time(1);
        for (;;) {
            // wait for 1 sec and check status
            auto status = future_.wait_for(wait_time);
            if (status == std::future_status::ready 
             || 0 == timeout_sec--)
                break;
        }
        return 0;
    }

    int write_dcr(uint32_t addr, uint32_t value) {
        if (future_.valid()) {
            future_.wait(); // ensure prior run completed
        }        
        processor_.write_dcr(addr, value);
        dcrs_.write(addr, value);
        return 0;
    }

    uint64_t read_dcr(uint32_t addr) const {
        return dcrs_.read(addr);
    }

private:

    RAM                 ram_;
    Processor           processor_;
    MemoryAllocator     global_mem_;
    MemoryAllocator     local_mem_;
    DeviceConfig        dcrs_;
    std::future<void>   future_;
};

///////////////////////////////////////////////////////////////////////////////

extern int vx_dev_caps(vx_device_h hdevice, uint32_t caps_id, uint64_t *value) {
   if (nullptr == hdevice)
        return  -1;

    vx_device *device = ((vx_device*)hdevice);

    switch (caps_id) {
    case VX_CAPS_VERSION:
        *value = IMPLEMENTATION_ID;
        break;
    case VX_CAPS_NUM_THREADS:
        *value = NUM_THREADS;
        break;
    case VX_CAPS_NUM_WARPS:
        *value = NUM_WARPS;
        break;
    case VX_CAPS_NUM_CORES:
        *value = NUM_CORES * NUM_CLUSTERS;
        break;
    case VX_CAPS_CACHE_LINE_SIZE:
        *value = CACHE_BLOCK_SIZE;
        break;
    case VX_CAPS_GLOBAL_MEM_SIZE:
        *value = GLOBAL_MEM_SIZE;
        break;
    case VX_CAPS_KERNEL_BASE_ADDR:
         *value = (uint64_t(device->read_dcr(VX_DCR_BASE_STARTUP_ADDR1)) << 32)
                          | device->read_dcr(VX_DCR_BASE_STARTUP_ADDR0);
        break;    
    case VX_CAPS_ISA_FLAGS:
        *value = ((uint64_t(MISA_EXT))<<32) | ((log2floor(XLEN)-4) << 30) | MISA_STD;
        break;
    default:
        std::cout << "invalid caps id: " << caps_id << std::endl;
        std::abort();
        return -1;
    }

    return 0;
}

extern int vx_dev_open(vx_device_h* hdevice) {
    if (nullptr == hdevice)
        return  -1;

    auto device = new vx_device();
    if (device == nullptr)
        return -1;

    int err = dcr_initialize(device);
    if (err != 0) {
        delete device;
        return err;
    }

#ifdef DUMP_PERF_STATS
    perf_add_device(device);
#endif

    *hdevice = device;

    return 0;
}

extern int vx_dev_close(vx_device_h hdevice) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    
#ifdef DUMP_PERF_STATS
    perf_remove_device(hdevice);
#endif

    delete device;

    return 0;
}

extern int vx_mem_alloc(vx_device_h hdevice, uint64_t size, int type, uint64_t* dev_addr) {
    if (nullptr == hdevice 
     || nullptr == dev_addr
     || 0 == size)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    return device->mem_alloc(size, type, dev_addr);
}

extern int vx_mem_free(vx_device_h hdevice, uint64_t dev_addr) {
    if (nullptr == hdevice)
        return -1;

    if (0 == dev_addr)
        return 0;

    vx_device *device = ((vx_device*)hdevice);
    return device->mem_free(dev_addr);
}

extern int vx_mem_info(vx_device_h hdevice, int type, uint64_t* mem_free, uint64_t* mem_used) {
    if (nullptr == hdevice)
        return -1;

    auto device = ((vx_device*)hdevice);
    return device->mem_info(type, mem_free, mem_used);
}

extern int vx_copy_to_dev(vx_device_h hdevice, uint64_t dev_addr, const void* host_ptr, uint64_t size) {
    if (nullptr == hdevice)
        return -1;

    auto device = (vx_device*)hdevice;
    return device->upload(dev_addr, host_ptr, size);
}

extern int vx_copy_from_dev(vx_device_h hdevice, void* host_ptr, uint64_t dev_addr, uint64_t size) {
    if (nullptr == hdevice)
        return -1;

    auto device = (vx_device*)hdevice;
    return device->download(host_ptr, dev_addr, size);
}

extern int vx_start(vx_device_h hdevice) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    return device->start();
}

extern int vx_ready_wait(vx_device_h hdevice, uint64_t timeout) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);
    return device->wait(timeout);
}

extern int vx_dcr_write(vx_device_h hdevice, uint32_t addr, uint64_t value) {
    if (nullptr == hdevice)
        return -1;

    vx_device *device = ((vx_device*)hdevice);

    // Ensure ready for new command
    if (vx_ready_wait(hdevice, -1) != 0)
        return -1;  
    return device->write_dcr(addr, value);
}