#ifndef LS_PRE_COMPILE_MACROS_HPP
#define LS_PRE_COMPILE_MACROS_HPP

#ifdef VIENNALS_USE_SHARED_LIBS

#define PRECOMPILE_PRECISION_DIMENSION(className)                              \
  typedef className<double, 2> className##_double_2;                           \
  typedef className<double, 3> className##_double_3;                           \
  typedef className<float, 2> className##_float_2;                             \
  typedef className<float, 3> className##_float_3;                             \
  extern template class className<double, 2>;                                  \
  extern template class className<double, 3>;                                  \
  extern template class className<float, 2>;                                   \
  extern template class className<float, 3>;

#else

// do nothing if we use header only
#define PRECOMPILE_PRECISION_DIMENSION(className)                              \
  typedef className<double, 2> className##_double_2;                           \
  typedef className<double, 3> className##_double_3;                           \
  typedef className<float, 2> className##_float_2;                             \
  typedef className<float, 3> className##_float_3;

#endif

#define PRECOMPILE_SPECIALIZE(className)                                       \
  template class className<double, 2>;                                         \
  template class className<double, 3>;                                         \
  template class className<float, 2>;                                          \
  template class className<float, 3>;

#endif // LS_PRE_COMPILE_MACROS_HPP
