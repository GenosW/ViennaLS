// STL
#include <iostream>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsTestAsserts.hpp>
//#define LS_TREE
#include <lsDimension.hpp>
#include <lsTest.hpp>
#include <lsVTKWriter.hpp>

template <int D = 3>
lsTestStatus testDimension(void)
{

  using Dim = Dimension<D>;
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

  return lsTestStatus::SUCCESS;
};

int main(int argc, char **argv)
{
  lsTest test_2d = INIT_LSTEST(test_2d);
  MAKE_LSTEST(test_3d);

  test_2d.run(testDimension<2>);
  test_3d.run(testDimension<3>);

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  check(test_2d, test_3d);

  if (all(test_2d.wasSuccess(), test_3d.wasSuccess()))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
}
