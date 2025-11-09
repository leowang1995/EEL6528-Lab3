
// save_detected_sample.cpp
// Saves detected blocks from the FIFO to samples_iq.bin
#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#include <signal.h>

#include "usrp_interface.hpp"

// These must be exposed from usrp_interface.cpp
extern BlockFIFO g_fifo;
extern std::atomic<bool> g_stop;

static void on_sigint(int) { 
    g_stop = true;
    std::cout << "\nStopping...\n";
}

int main() {
    signal(SIGINT, on_sigint);
    
    std::ofstream ofs("samples_iq.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Cannot open samples_iq.bin\n";
        return 1;
    }
    
    std::cout << "Saving to samples_iq.bin (Ctrl+C to stop)...\n";
    
    std::vector<std::complex<float>> blk;
    int count = 0;
    
    while (!g_stop) {
        if (!g_fifo.pop(blk)) continue;
        
        // Write block to file
        ofs.write(reinterpret_cast<const char*>(blk.data()),
                  blk.size() * sizeof(std::complex<float>));
        
        if (++count % 10 == 0) {
            std::cout << "Blocks saved: " << count << "\r" << std::flush;
        }
    }
    
    ofs.close();
    std::cout << "\nTotal blocks: " << count << "\n";
    return 0;
}
