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

                std::function<void()> task;
                { //New scope: Lock and wait for action
                    std::cout << "Waiting" << std::endl;
                    std::unique_lock<std::mutex> lock(mtx_queue);
                    event_cond.wait(lock, [=] { return !tasks.empty() || request_stop; });
                    if (request_stop)
                        break;

                    task = std::move(tasks.front());
                    tasks.pop();

                }
                task();
            }



        });
    }
}

void thread_pool::stop() {
    {
        std::unique_lock<std::mutex> lock(mtx_queue);
        request_stop = true;
    }

    event_cond.notify_all();
    for (auto &t : threads) {
        t.join();
    }

}



void thread_pool::enqueue(std::function<void()> task) {
    {
        std::unique_lock<std::mutex> lock(mtx_queue);
        tasks.emplace(std::move(task));
    }
    event_cond.notify_one();
}
