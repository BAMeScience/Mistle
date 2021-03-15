#ifndef SIMPLE_EXAMPLE_THREAD_POOL_H
#define SIMPLE_EXAMPLE_THREAD_POOL_H

#include <functional>
#include <thread>
#include <vector>
#include <condition_variable>
#include <cstddef>

class thread_pool {

    std::vector<std::thread> threads;
    //todo add task vector std::vector<std::future>
    std::mutex mtx;
    std::condition_variable event_cond;
    bool request_stop = false;

    size_t size;



public:
    thread_pool(size_t n);
    ~thread_pool();
    void start();
    void stop();
    void wait();
};


#endif //SIMPLE_EXAMPLE_THREAD_POOL_H
