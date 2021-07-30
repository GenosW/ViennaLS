// STL
#include <iostream>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
#include <lsTestAsserts.hpp>
//#define LS_TREE
#include <lsDimension.hpp>
#include <lsTree.hpp>
#include <lsTest.hpp>
#include <lsVTKWriter.hpp>

template <int D>
void checkTree(lsTree<double, D> treeToCheck, size_t N)
{
  auto numLevels = treeToCheck.getDepth();
  using Dim = Dimension<D>;
  Dim dim = Dim(D); // start with highest dim
  std::vector<Dim> orderOfDims(numLevels + 1, 0);
  for (auto &dimOfLevel : orderOfDims)
  {
    dimOfLevel = dim;
    // std::cout << dimOfLevel << " -> ";
    --dim;
  }

  std::cout << std::endl;

  std::vector<size_t> ptsInTree(numLevels + 1, 0);
  for (auto &node : treeToCheck.getTreeNodes())
  {
    // std::cout << node->identifier << node->dimSplit << " ?= " << orderOfDims[node->level] << std::endl;
    // std::cout << node->dimSplit << " ?= " << orderOfDims[node->level] << std::endl;
    ptsInTree[node->level] += node->size();
    // Check if node was split in correct dimension
    LSTEST_ASSERT(node->dimSplit == orderOfDims[node->level]);
  }
  // Check if all points were sorted and are present in each level of tree
  LSTEST_ASSERT(
      std::all_of(ptsInTree.cbegin(), ptsInTree.cend(),
                  [N](auto item)
                  { return item == N; }));
};

template <int D>
lsTestStatus testTree(void)
{
  //constexpr int D = Dim;
  std::cout << "D: " << D << std::endl;
  const std::string prefix = "meshes/" + std::to_string(D) + "D";

  setup_omp(4);

  // Setup reference geometry
  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const double radius = 10;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();
  size_t N = levelSet->getDomain().getNumberOfPoints();

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << N << std::endl;

  // --- Volumetric mesh
  std::cout << "--- Mesh: " << std::endl;
  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, prefix + "_Mesh.vtk").apply();

  // Get reference geometry parameters
  LSTEST_ASSERT(mesh->getNodes().size() == N);

  // Compute lsTree of volumetric mesh
  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh.vtu").apply();

  // tree.printBT();
  // EVALUATE

  checkTree(tree, N);

  //---DiskMesh
  std::cout << "--- DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_DiskMesh.vtu").apply();

  size_t N2 = mesh->getNodes().size();
  std::cout << "N2: " << N2 << std::endl;

  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh.vtu").apply();

  checkTree(tree2, N2);

  // PASS!
  return lsTestStatus::SUCCESS;
};

int main(int argc, char **argv)
{
  lsTest test_2d = INIT_LSTEST(test_2d);
  MAKE_LSTEST(test_3d);

  test_2d.run(testTree<2>);
  std::cout << "------------- NEXT TEST -------------" << std::endl;
  test_3d.run(testTree<3>);

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  check(test_2d, test_3d);

  if (all(test_2d.wasSuccess(), test_3d.wasSuccess()))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
}
