#ifndef LS_TEST
#define LS_TEST
#include <iostream>
#define PRINT(Var) std::cout << #Var << ": " << Var

enum struct lsTestStatus : unsigned
{
    SUCCESS = 0,
    FAILED = 1,
    UNCHECKED = 2

};

struct lsTest
{
    lsTestStatus status = lsTestStatus::UNCHECKED;
    std::string name = "unnamed";

    lsTest(const lsTestStatus stat)
    {
        status = stat;
    }
    lsTest(lsTestStatus stat, std::string nam) : status(stat), name(nam) {}

    template <typename T>
    bool run(T (&test)(void))
    {
        this->status = lsTestStatus{(uint)test()};
        return this->status == lsTestStatus::SUCCESS;
    }

    void check() const
    {
        std::cout << name;
        if (status == lsTestStatus::FAILED)
        {
            std::cout << ": FAILED";
        }
        else if (status == lsTestStatus::SUCCESS)
        {
            std::cout << ": SUCCESS";
        }
        else if (status == lsTestStatus::UNCHECKED)
        {
            std::cout << ": UNCHECKED";
        }
        else
        {
            std::cout << ": unknown status";
        }
        std::cout << std::endl;
    };

    lsTestStatus operator()() const
    {
        return this->status;
    }

    bool wasSuccess() const
    {
        return this->status == lsTestStatus::SUCCESS;
    }
};

#define INIT_LSTEST(Var) lsTest(lsTestStatus::UNCHECKED, #Var)
#define MAKE_LSTEST(Var) lsTest Var(lsTestStatus::UNCHECKED, #Var)

template <class Test, class... Tests>
void check(Test const &t1, Tests const &...rest)
{
    t1.check();
    if constexpr (sizeof...(rest) > 0)
    {
        check(rest...);
    }
}

template <typename... Args>
bool all(Args... args)
{
    return (... && args);
}

template <typename T, typename... Args>
bool all_equal(T target, Args... args)
{
    return ((target == args) && ...);
}

void setup_omp(uint numThreads)
{
    int maxThreads = omp_get_max_threads();
    omp_set_num_threads(4);
    std::cout << "Using " << omp_get_max_threads() << " (of max " << maxThreads << ") threads for this test." << std::endl;
};
#endif // LS_TEST