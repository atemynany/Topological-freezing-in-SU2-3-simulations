// ==============================================================================
// progressbar.hh
// ==============================================================================
// Simple console progress bar for tracking simulation progress.
//
// Author: Alexander de Barros Noll
// Date: January 2026
// ==============================================================================

#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH

#include <iostream>

/**
 * @brief Display a progress bar in the console.
 * 
 * @param progress Progress value between 0.0 and 1.0
 * @param bar_width Width of the progress bar in characters (default: 50)
 */
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

/**
 * @brief Clear the progress bar line.
 * 
 * @param bar_width Width of the progress bar in characters (default: 50)
 */
inline void progress_bar_clear(int bar_width = 50) {
    std::cout << "\r";
    for (int i = 0; i < bar_width + 10; ++i)
        std::cout << " ";
    std::cout << "\r";
    std::cout.flush();
}

#endif // PROGRESSBAR_HH
