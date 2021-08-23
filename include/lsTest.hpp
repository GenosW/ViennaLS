#ifndef LS_TEST
#define LS_TEST
#include <iostream>
#define PRINT(Var) std::cout << #Var << ": " << Var

namespace lsInternal
{
  template <typename... Args>
  bool all(Args... args) { return (... && args); }

  template <typename T, typename... Args>
  bool all_equal(T target, Args... args)
  {
    return ((target == args) && ...);
  }

  void setup_omp(uint numThreads)
  {
    int maxThreads = omp_get_max_threads();
    omp_set_num_threads(4);
    std::cout << "Using " << omp_get_max_threads() << " (of max " << maxThreads
              << ") threads for this test." << std::endl;
  };
}
#endif // LS_TEST