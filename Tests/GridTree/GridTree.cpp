// STL
#include <iostream>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsDimension.hpp>
#include <lsTest.hpp>
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>

template <int D>
uint checkTree(lsTree<double, D> treeToCheck, size_t N)
{
  auto numLevels = treeToCheck.getDepth();
  using Dim = typename lsInternal::Dimension<D>;
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
    // std::cout << node->identifier << node->dimensionToSplit << " ?= " <<
    // orderOfDims[node->level] << std::endl; std::cout <<
    // node->dimensionToSplit << " ?= " << orderOfDims[node->level] <<
    // std::endl;
    ptsInTree[node->level] += node->size();
    // Check if node was split in correct dimension
    LSTEST_ASSERT(node->dimensionToSplit == orderOfDims[node->level]);
  }
  // Check if all points were sorted and are present in each level of tree
  LSTEST_ASSERT(std::all_of(ptsInTree.cbegin(), ptsInTree.cend(),
                            [N](auto item)
                            { return item == N; }));
  return 1;
};

template <class T, int D, class container>
void queryTree(lsTree<T, D> tree, container point)
{

  auto neighborhood = tree.getBin(point);
  std::cout << "leafNum: " << neighborhood->leafNum << std::endl;
  std::cout << "color: " << neighborhood->color << std::endl;
  std::cout << "level: " << neighborhood->level << std::endl;
  std::cout << "identifier: " << neighborhood->identifier << std::endl;
  std::cout << "median: " << neighborhood->dimensionToSplit << " = " << neighborhood->median << std::endl;
}

template <int D>
uint testTree(void)
{
  // constexpr int D = Dim;
  std::cout << "D: " << D << std::endl;
  const std::string prefix = "meshes/sphere/" + std::to_string(D) + "D";

  lsInternal::setup_omp(4);

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
  lsVTKWriter<double>(mesh, lsFileFormatEnum::VTU, prefix + "_Mesh").apply();

  // Get reference geometry parameters
  LSTEST_ASSERT(mesh->getNodes().size() == N);

  // Compute lsTree of volumetric mesh
  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  tree.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh").apply();

  // tree.printBT();
  // EVALUATE

  checkTree(tree, N);

  //---DiskMesh
  std::cout << "--- DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_DiskMesh").apply();

  size_t N2 = mesh->getNodes().size();
  std::cout << "N2: " << N2 << std::endl;

  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  tree2.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh").apply();

  checkTree(tree2, N2);

  auto ptForQuery = centre;
  ptForQuery[0] += radius;
  queryTree(tree, ptForQuery);

  // PASS!
  return 1;
};

template <int D>
uint testTreeWithBox(void)
{
  std::cout << "D: " << D << std::endl;
  const std::string prefix = "meshes/box/" + std::to_string(D) + "D";

  lsInternal::setup_omp(4);

  // Setup reference geometry
  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  double a = 30, b = 20, c = 10;
  double bounds[2 * 3] = {-a, a, -b, b, 0, 0};
  if constexpr (D == 3)
  {
    bounds[4] = -c;
    bounds[5] = c;
  }

  hrleVectorType<double, D> botLeft(bounds[0], bounds[2], bounds[4]);
  hrleVectorType<double, D> topRight(bounds[1], bounds[3], bounds[5]);
  lsMakeGeometry<double, D>(
      levelSet, lsSmartPointer<lsBox<double, D>>::New(botLeft, topRight))
      .apply();
  size_t N = levelSet->getDomain().getNumberOfPoints();

  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << N << std::endl;

  // --- Volumetric mesh
  std::cout << "--- Mesh: " << std::endl;
  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, lsFileFormatEnum::VTU, prefix + "_Mesh").apply();

  // Get reference geometry parameters
  LSTEST_ASSERT(mesh->getNodes().size() == N);

  // Compute lsTree of volumetric mesh
  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  tree.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh").apply();

  // tree.printBT();
  // EVALUATE

  checkTree(tree, N);

  //---DiskMesh
  std::cout << "--- DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_DiskMesh").apply();

  size_t N2 = mesh->getNodes().size();
  std::cout << "N2: " << N2 << std::endl;

  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  tree2.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh").apply();

  checkTree(tree2, N2);

  queryTree(tree, botLeft);

  // PASS!
  return 1;
};

int main(int argc, char **argv)
{
  uint test_2d_sphere = 0;
  uint test_3d_sphere = 0;
  uint test_2d_box = 0;
  uint test_3d_box = 0;

  std::cout << "############### SPHERE ##############" << std::endl;
  test_2d_sphere = testTree<2>();
  std::cout << "------------- NEXT TEST -------------" << std::endl;
  test_3d_sphere = testTree<3>();
  std::cout << "################ BOX ################" << std::endl;
  test_2d_box = testTreeWithBox<2>();
  std::cout << "------------- NEXT TEST -------------" << std::endl;
  test_3d_box = testTreeWithBox<3>();

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  if (lsInternal::all(test_2d_sphere, test_3d_sphere,
                      test_2d_box, test_3d_box))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
};
