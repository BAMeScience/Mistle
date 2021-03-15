#ifndef SIMPLE_EXAMPLE_THREAD_POOL_H
#define SIMPLE_EXAMPLE_THREAD_POOL_H

#include <functional>
#include <thread>
#include <vector>
#include <condition_variable>
#include <cstddef>

class thread_pool {

    std::vector<std::thread> threads;

    thread_pool(size_t n);

};


#endif //SIMPLE_EXAMPLE_THREAD_POOL_H
