// STL
#include <iostream>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsTestAsserts.hpp>
#include <lsDimension.hpp>
#include <lsTest.hpp>
#include <lsVTKWriter.hpp>

template <int D = 3>
uint testDimension(void)
{

  using Dim = typename lsInternal::Dimension<D>;
  std::vector<uint> shouldBe;

  // test forward iteration
  Dim dim = Dim(1);
  std::cout << "Iterate (++): " << dim;
  if constexpr (D == 3)
    shouldBe = std::vector<uint>{1, 2, 0};
  if constexpr (D == 2)
    shouldBe = std::vector<uint>{1, 0};
  for (auto i : shouldBe)
  {
    std::cout << " --> " << ++dim << std::flush;
    uint arrIdx = dim.toArrayIndex();
    LSTEST_ASSERT(arrIdx == i);
  }
  std::cout << std::endl;

  // test forward copy iteration
  Dim dim2 = dim++;
  std::cout << "Copy (+): " << dim << " vs " << dim2 << std::endl;
  LSTEST_ASSERT(dim2 != dim);

  // test backward iteration
  dim = Dim(D);
  std::cout << "Iterate (--): " << dim;
  if constexpr (D == 3)
    shouldBe = std::vector<uint>{1, 0, 2};
  if constexpr (D == 2)
    shouldBe = std::vector<uint>{0, 1};
  for (auto i : shouldBe)
  {
    std::cout << " --> " << --dim << std::flush;
    uint arrIdx = dim.toArrayIndex();
    LSTEST_ASSERT(arrIdx == i);
  }
  std::cout << std::endl;

  // test backward copy iteration
  dim2 = dim--;
  std::cout << "Copy (-): " << dim << " vs " << dim2 << std::endl;
  LSTEST_ASSERT(dim2 != dim);

  uint initDim = 2;
  std::cout << "type Conversions" << std::endl;
  std::cout << "Assign from uint (expected: y == 1)" << std::endl;
  Dim dim3 = initDim; // y

  uint shouldBeThis = 1; // minus 1 included
  std::cout << "implicit to uint: " << dim3 << " == " << shouldBeThis
            << std::endl;
  LSTEST_ASSERT(dim3 == shouldBeThis)
  std::cout << "explicit to uint: " << dim3 << " == " << shouldBeThis
            << std::endl;
  LSTEST_ASSERT(uint(dim3) == shouldBeThis)
  std::cout << "explicit to int: " << dim3 << " == " << shouldBeThis
            << std::endl;
  LSTEST_ASSERT(int(dim3) == shouldBeThis)
  std::cout << "explicit to size_t: " << dim3 << " == " << shouldBeThis
            << std::endl;
  LSTEST_ASSERT(size_t(dim3) == shouldBeThis)

  return 1;
};

int main(int argc, char **argv)
{
  uint test_2d, test_3d;

  test_2d = testDimension<2>();
  test_3d = testDimension<3>();

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  if (lsInternal::all(test_2d, test_3d))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
}
