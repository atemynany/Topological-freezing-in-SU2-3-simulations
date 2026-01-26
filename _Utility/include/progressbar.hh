// progressbar.hh
#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH

#include <iostream>

inline void progress_bar(double progress, int bar_width = 50) {
    std::cout << "[";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

inline void progress_bar_clear(int bar_width = 50) {
    std::cout << "\r";
    for (int i = 0; i < bar_width + 10; ++i)
        std::cout << " ";
    std::cout << "\r";
    std::cout.flush();
}

#endif