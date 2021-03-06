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

    std::mutex mtx_queue;
    std::condition_variable event_cond;
    std::condition_variable finished_cond;
    bool request_stop = false;

    size_t size;
    int busy_threads = 0;

public:
    explicit thread_pool(size_t n);
    ~thread_pool();
    void start();
    void stop();


    void add_thread();
    void wait_for_all_threads();
    void join_all();

    //template<typename F, typename... Args>
    //void enqueue(F f, Args&&... args);

    void enqueue(std::function<void()> task);
    std::mutex mtx;


    void thread_waiting_loop();

    //Getter
    size_t get_size() const;
private:

};


#endif //SIMPLE_EXAMPLE_THREAD_POOL_H
