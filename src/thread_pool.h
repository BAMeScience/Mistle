#ifndef SIMPLE_EXAMPLE_THREAD_POOL_H
#define SIMPLE_EXAMPLE_THREAD_POOL_H

#include <functional>
#include <future>
#include <thread>
#include <vector>
#include <queue>
#include <condition_variable>
#include <cstddef>

class thread_pool {



    std::vector<std::thread> threads;
    std::queue<std::function<void()>> tasks;
    //todo add task vector std::vector<std::future>
    std::mutex mtx_queue;
    std::condition_variable event_cond;
    bool request_stop = false;

    size_t size;



public:
    thread_pool(size_t n);
    ~thread_pool();
    void start();
    void stop();
    void wait();

    void enqueue(std::function<void()> task);

    std::mutex mtx;

};


#endif //SIMPLE_EXAMPLE_THREAD_POOL_H
