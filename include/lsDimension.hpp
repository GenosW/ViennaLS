#ifndef LS_DIMENSION
#define LS_DIMENSION

//#include <iostream>
#include <ostream>

namespace lsInternal
{
  template <int D = 3>
  struct Dimension
  {
    using enum_type = uint;
    enum options : enum_type
    {
      START = 0,
      x = 1,
      y = 2,
      z,
      END = D + 1
    };
    options value = options(1);

    Dimension() = default;
    Dimension(Dimension const &) = default;
    Dimension(enum_type init) : value(options(init + 1)) {}

    Dimension &operator=(Dimension other)
    {
      this->value = other.value;
      return *this;
    }

    template <typename T = uint>
    T operator()() const
    {
      return static_cast<T>(value);
    };

    // std::ostream &operator<<(std::ostream &os, const Dimension<D> &dim);

    Dimension &operator++()
    {
      value = static_cast<options>(static_cast<options>(value) + 1);
      if (value == END)
        value = x;
      return *this;
    }

    Dimension operator++(int)
    {
      Dimension new_dim = *this;
      ++new_dim;
      return new_dim;
    }

    Dimension &operator--()
    {
      value = static_cast<options>(static_cast<options>(value) - 1);
      if (value == START)
        value = options(D);
      return *this;
    }

    Dimension operator--(int)
    {
      Dimension new_dim = *this;
      --new_dim;
      return new_dim;
    }

    bool operator==(Dimension &other) { return other.value == this->value; }

    bool operator!=(Dimension &other) { return other.value != this->value; }

    operator enum_type() const { return toArrayIndex<enum_type>(); }

    template <typename T = size_t>
    T toArrayIndex() const
    {
      return static_cast<T>(value) - 1;
    }
  };

  template <int D>
  std::ostream &operator<<(std::ostream &os, const Dimension<D> &dim)
  {
    switch (dim())
    {
    case Dimension<D>::x:
      os << "x";
      break;
    case Dimension<D>::y:
      os << "y";
      break;
    case Dimension<D>::z:
      os << "z";
      break;
    default:
      os.setstate(std::ios_base::failbit);
    }
    return os;
  };
} // end namespace lsInternal

#endif // LS_DIMENSION