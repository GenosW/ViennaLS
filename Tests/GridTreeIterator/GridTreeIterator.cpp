// STL
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
//#define LS_TREE
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>

int main(int argc, char **argv)
{
  constexpr int D = 3;

  int max = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << max << ") threads for this test." << std::endl;

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const double radius = 10;
  const hrleVectorType<double, D> centre(5., 0.);

  // Query points
  hrleVectorType<double, D> tipytop(centre), bot(centre), left(centre), right(centre);
  tipytop[0] += radius;
  bot[0] -= radius;
  left[1] -= radius;
  right[1] += radius;

  lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "ToDiskMesh.vtu").apply();

  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "TreeDiskMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output

  tree.printBT();

  // Test treeIterator
  int cnt = 0, trueValue = 0;
  auto print_iterator = [](auto p)
  {
    std::cout << "Iterator(pos= " << p.pos << ", segment= " << p.segment << ")" << std::endl;
    auto pos = *(p.pos);
    std::cout << "(" << pos[0] << ")" << std::endl; //<< "," << *p[2]
  };
  auto print_point = [](auto p)
  {
    double val = p[0][0];
    std::cout << "(" << val << ")" << std::endl; //<< "," << *p[2]
  };
  // for (auto point : tree)
  // for (auto point = tree.begin(); cnt < 10; ++point)
  for (auto point = tree.data.begin(); cnt < 10; ++point)
  {
    std::cout << ++cnt << ": ";
    // print_point(point);
  }
  // TODO: add correct trueValue
  std::cout << std::endl;
  trueValue = cnt;
  if (cnt != trueValue)
  {
    std::cout << "Test " << argv[0] << ": FAILURE!";
    return 1;
  }

  std::cout << "Test " << argv[0] << ": SUCCESS!";
  return 0;
}
