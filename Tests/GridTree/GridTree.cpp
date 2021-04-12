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

int main(int argc, char** argv) {

  constexpr int D = 2;

  int before = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << before << ") threads for this test." << std::endl;

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();
  
  const double radius = 10;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, 2>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "ToMesh.vtk").apply();

  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "TreeMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
  tree.printInfo();
  tree.printTree();
  tree.printTreeByLevel();

  std::cout << "DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "ToDiskMesh.vtu").apply();

  std::cout << "VTKWriter --> Disk.vtu: DONE" << std::endl;
  tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "TreeDiskMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
  tree.printInfo();
  tree.printTree();
  tree.printTreeByLevel();
  
  std::cout << "Test " << argv[0] << ": SUCCESS!";
  return 0;
}
