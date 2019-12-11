/**
 *  @file    ExecutionTimer.hpp
 *  @brief   Measure the elapsed time of a code block.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef EXECUTION_TIMER_H
#define EXECUTION_TIMER_H
#include <chrono>
#include "spdlog/spdlog.h"

namespace utils {

/**
* Measure the elapsed time of a code block. The default constructor starts
* the timer.
* @param T The units -- default ms.
*/
template<typename T = std::chrono::milliseconds>
class ExecutionTimer {
private:
    /**
    * The reference clock count.
    */
    const std::chrono::steady_clock::time_point m_start = 
      std::chrono::steady_clock::now();

public:
    ExecutionTimer() = default;
    ~ExecutionTimer() {
        const auto end = std::chrono::steady_clock::now();
        int out = std::chrono::duration_cast<T>( end - m_start ).count();
        spdlog::debug("Destructor Elapsed: {:08d}", out);
    }
    /**
    * Returns the thick count.
    * @return the tick count in selected units.
    */
    inline unsigned int stop() {
        const auto end = std::chrono::steady_clock::now();
        unsigned int out = 
            std::chrono::duration_cast<T>( end - m_start ).count();
        spdlog::debug("Stop Elapsed: {:8d}", out);
        return out;
    }

}; // ExecutionTimer

}  // utils

#endif // EXECUTION_TIMER_H