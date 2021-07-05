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

int test_2D(void)
{

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
  lsVTKWriter<double>(mesh, "2D_ToMesh.vtk").apply();

  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "2D_TreeMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
  tree.printInfo();
  // tree.printTree();
  // tree.printTreeByLevel();
  tree.printBT();

  std::cout << "DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "2D_ToDiskMesh.vtu").apply();

  std::cout << "VTKWriter --> Disk.vtu: DONE" << std::endl;
  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "2D_TreeDiskMesh.vtu").apply();
  // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
  tree2.printInfo();
  // tree2.printTree();
  // tree2.printTreeByLevel();
  tree.printBT();
  return EXIT_SUCCESS;
};

int test_3D(void)
{
  constexpr int D = 3;

  int before = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << before << ") threads for this test." << std::endl;

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const double radius = 3;
  for (double radius = 1; radius < 6.1; ++radius)
  {
    // radius = 1 --> 25 po{ints
    // radius = 2 --> 86 points‚
    // radius = 3 --> 230 points
    // radius = 4 -->  points
    const hrleVectorType<double, D> centre(5., 0., 0.);

    lsMakeGeometry<double, D>(
        levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
        .apply();

    std::cout << "Initial: " << std::endl;
    std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
              << std::endl;
    lsToMesh<double, D>(levelSet, mesh).apply();
    lsVTKWriter<double>(mesh, "3D_ToMesh.vtk").apply();
#include <chrono>
    const auto start = std::chrono::high_resolution_clock::now();
    auto tree = lsTree<double, D>(mesh);
    tree.apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTU, "3D_TreeMesh.vtu").apply();
    // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
    tree.printInfo();
    // tree.printTree();
    // tree.printTreeByLevel();
    tree.printBT();

    std::cout << "DiskMesh: " << std::endl;
    lsToDiskMesh<double, D>(levelSet, mesh).apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTU, "3D_ToDiskMesh.vtu").apply();

    std::cout << "VTKWriter --> Disk.vtu: DONE" << std::endl;
    auto tree2 = lsTree<double, D>(mesh);
    tree2.apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTU, "3D_TreeDiskMesh.vtu").apply();
    // TODO: remove verbose test output and replace with SUCCESS/FAIL style output
    tree2.printInfo();
    // tree2.printTree();
    // tree2.printTreeByLevel();
    tree.printBT();
  }
  return EXIT_SUCCESS;
};

int main(int argc, char **argv)
{
  int success_2d_test = EXIT_SUCCESS, success_3d_test = EXIT_SUCCESS;
  // success_2d_test = test_2D();
  success_3d_test = test_3D();
  if (success_2d_test == EXIT_SUCCESS && success_3d_test == EXIT_SUCCESS)
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return 0;
}
