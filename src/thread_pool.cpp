#include <iostream>
#include "thread_pool.h"

thread_pool::thread_pool(size_t n) : size(n) {
    start();
}

thread_pool::~thread_pool() {
    stop();
}

void thread_pool::wait() {

}

void thread_pool::start() {
    for (size_t i = 0; i < size; ++i) {
        threads.emplace_back([=] {
            while (true) {
                std::cout << "Waiting" << std::endl;
                std::unique_lock<std::mutex> lock(mtx);
                event_cond.wait(lock, [=] {return request_stop;});
                if (request_stop)
                    break;
            }
        });
    }
}

void thread_pool::stop() {
    {
        std::unique_lock<std::mutex> lock(mtx);
        request_stop = true;
    }

    event_cond.notify_all();
    for (auto &t : threads) {
        t.join();
    }

}
