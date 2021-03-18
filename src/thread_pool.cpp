#include <iostream>
#include "thread_pool.h"

thread_pool::thread_pool(size_t n) : size(n) {
    start();
    busy_threads = 0;
}

thread_pool::~thread_pool() {
    stop();
}

void thread_pool::start() {
    for (size_t i = 0; i < size; ++i) {
        threads.emplace_back(&thread_pool::thread_waiting_loop, this); //pushback std::thread
    }
}

void thread_pool::stop() {
    {
        std::unique_lock<std::mutex> lock(mtx_queue);
        request_stop = true;
    }

    event_cond.notify_all();
    join_all();

}



void thread_pool::enqueue(std::function<void()> task) {
    {
        std::unique_lock<std::mutex> lock(mtx_queue);
        tasks.emplace(task);
    }
    event_cond.notify_one();
}

void thread_pool::join_all() {
    for (auto &t : threads) {
        t.join();
    }
}

void thread_pool::wait_for_all_threads() {
    std::unique_lock<std::mutex> lock(mtx_queue);
    finished_cond.wait(lock, [this] { return (tasks.empty() && (busy_threads == 0)); });
}

void thread_pool::thread_waiting_loop() {
    while (true) {

        std::function<void()> task;
        { //New scope: Lock and wait for action

            std::unique_lock<std::mutex> lock(mtx_queue); //TODO check std:acquire defer lock
            event_cond.wait(lock, [this] { return !tasks.empty() || request_stop; });
            if (request_stop && tasks.empty())
                break;

            ++busy_threads;
            task = std::move(tasks.front());
            tasks.pop();

        }
        task();

        mtx_queue.lock();
        --busy_threads;
        finished_cond.notify_one();
        mtx_queue.unlock();

    }
}

/*
template<typename F, typename... Args>
void thread_pool::enqueue(F f, Args &&... args) {
    tasks.emplace( std::bind(f, std::forward<Args>(args)...) );
}*/
