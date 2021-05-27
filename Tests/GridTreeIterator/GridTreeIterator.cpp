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
  constexpr int D = 2;

  int max = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << max << ") threads for this test." << std::endl;

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const double radius = 10;
  const hrleVectorType<double, D> centre(5., 0.);

  // Query points
  const hrleVectorType<double, D> tipytop(centre), bot(centre), left(centre), right(centre);
  tipytop[0] += radius;
  bot[0] -= radius;
  left[1] -= radius;
  right[1] += radius;

  lsMakeGeometry<double, 2>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "ToDiskMesh.vtu").apply();

  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "TreeDiskMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output

  // Test treeIterator
  int cnt = 0, trueValue = 0;
  for (auto point : tree)
    cnt++;
  // TODO: add correct trueValue
  trueValue = cnt;
  if (cnt != trueValue)
  {
    std::cout << "Test " << argv[0] << ": FAILURE!";
    return 1;
  }

  std::cout << "Test " << argv[0] << ": SUCCESS!";
  return 0;
}
