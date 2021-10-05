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

#define DEBUG_INFO

template <typename T>
bool isInBounds(const T &value, const T &low, const T &high)
{
  return (value >= low) && (value < high);
}

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

template <class T, int D>
void showTreeRange(lsTree<T, D> &tree)
{
  auto &treeNodes = tree.getTreeNodes();
  auto size = treeNodes.size();
  auto numberOfLeafs = tree.getNumberOfLeafs();
  auto &first = treeNodes[size - numberOfLeafs];
  auto &last = treeNodes.back();
  LSTEST_ASSERT(first->isLeaf);
  std::cout << "numberOfLeafs:   [" << numberOfLeafs << "]\n";
  std::cout << "colorRange:   [" << first->color << ", " << last->color << "]\n";
  std::cout << "nodeIdxRange: [" << first->leafNum << ", " << last->leafNum << "]" << std::endl;
}

template <class nodePointer>
void printNodeInfo(nodePointer node)
{
  std::cout << "    leafNum: " << node->leafNum << "\n";
  std::cout << "    color: " << node->color << "\n";
  std::cout << "    level: " << node->level << "\n";
  std::cout << "    identifier: " << node->identifier << "\n";
  std::cout << "    median: " << node->dimensionToSplit << " = " << node->median << std::endl;
}

template <class T, int D, class container>
void queryTree(lsTree<T, D> tree, container point, int expectedBot, int expectedTop, std::string name)
{

  bool trace = true;
#ifdef DEBUG_INFO
  trace = true;
#endif
  std::cout << "Getting neighborhood for point: " << name << point << "\n";
  auto neighborhood = tree.getBin(point, trace);
  printNodeInfo(neighborhood);
  if (expectedBot != 0 && expectedTop != 0)
  {
    bool checkResult = isInBounds(neighborhood->color, expectedBot, expectedTop);
    std::cout << "    ("
              << expectedBot << " <= " << neighborhood->color
              << " <= " << expectedTop << "): " << checkResult << std::endl;
  }
}

template <int D>
uint testTreeWithSphere(void)
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
  // lsVTKWriter<double>(mesh, lsFileFormatEnum::VTU, prefix + "_Mesh").apply();

  // Get reference geometry parameters
  LSTEST_ASSERT(mesh->getNodes().size() == N);

  auto left = centre, top = centre;
  left[0] -= radius;
  top[D - 1] += radius;

  // Compute lsTree of volumetric mesh
  //   lsTree<double, D> tree(mesh);
  //   // tree.setMeanMethod(lsInternal::meanEnum::HARMONIC);
  //   tree.apply();
  //   // Debug tools
  //   tree.addColor(mesh); // TODO: Wrap in creation of mesh based on lsTree.data nodes
  //   lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh").apply();
  // #ifdef DEBUG_INFO
  //   tree.printBT();
  //   tree.dumpBins();
  // #endif
  //   // EVALUATE

  //   checkTree(tree, N);

  //   showTreeRange(tree);

  //   queryTree(tree, top, 7, 7, "top");
  //   queryTree(tree, left, -11, -11, "left");

  //   using treeNode = typename lsTree<double, D>::treeNode;

  //   auto print_node = [&tree](lsSmartPointer<treeNode> &bin)
  //   {
  //     for (auto it = bin->begin(tree); it != bin->end(tree); ++it)
  //     {
  //       lsInternal::print_point(*it);
  //       std::cout << endl;
  //     };
  //   };

  // lsSmartPointer<treeNode> bin = tree.getBin(top);
  // print_node(bin);

  // auto &treeNodes = tree.getTreeNodes();
  // auto size = treeNodes.size();
  // auto numberOfLeafs = tree.getNumberOfLeafs();

  // for (auto binIterator = treeNodes.begin() + (size - numberOfLeafs); binIterator != treeNodes.end(); ++binIterator)
  // {
  //   lsSmartPointer<treeNode> bin = lsSmartPointer<treeNode>(*binIterator);
  //   std::cout << "bin" << bin->leafNum << std::endl;
  //   print_node(bin);
  //   printNodeInfo(bin);
  // }

  // //---DiskMesh
  std::cout << "--- DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();

  size_t N2 = mesh->getNodes().size();
  std::cout << "N2: " << N2 << std::endl;

  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  tree2.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh").apply();
#ifdef DEBUG_INFO
  tree2.printBT();
  tree2.dumpBins();
#endif

  checkTree(tree2, N2);

  showTreeRange(tree2);
  queryTree(tree2, top, -15, 15, "top");
  queryTree(tree2, left, -15, -15, "left");

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

  std::cout << "Box: \n";
  std::cout << "         #--------------------------" << topRight << " :topRight"
            << "\n";
  std::cout << "botLeft: " << botLeft << "-----------------------#"
            << "\n";
  std::cout << "N: " << N << std::endl;

  // --- Volumetric mesh
  std::cout << "--- Mesh: " << std::endl;
  lsToMesh<double, D>(levelSet, mesh).apply();
  // lsVTKWriter<double>(mesh, lsFileFormatEnum::VTU, prefix + "_Mesh").apply();

  // Get reference geometry parameters
  LSTEST_ASSERT(mesh->getNodes().size() == N);

  // Compute lsTree of volumetric mesh
  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  tree.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh").apply();
#ifdef DEBUG_INFO
  tree.printBT();
  tree.dumpBins();
#endif

  // EVALUATE

  checkTree(tree, N);

  showTreeRange(tree);
  queryTree(tree, botLeft, -15, 15, "botLeft");
  queryTree(tree, topRight, -15, 15, "topRight");

  //---DiskMesh
  std::cout << "--- DiskMesh: " << std::endl;
  lsToDiskMesh<double, D>(levelSet, mesh).apply();
  // lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_DiskMesh").apply();

  size_t N2 = mesh->getNodes().size();
  std::cout << "N2: " << N2 << std::endl;

  auto tree2 = lsTree<double, D>(mesh);
  tree2.apply();
  tree2.addColor(mesh);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh").apply();
#ifdef DEBUG_INFO
  tree2.printBT();
  tree2.dumpBins();
#endif

  checkTree(tree2, N2);

  showTreeRange(tree2);
  queryTree(tree2, botLeft, -15, 15, "botLeft");
  queryTree(tree2, topRight, -15, 15, "topRight");

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
  test_2d_sphere = testTreeWithSphere<2>();
  std::cout << "------------- NEXT TEST -------------" << std::endl;
  // test_3d_sphere = testTreeWithSphere<3>();
  // std::cout << "################ BOX ################" << std::endl;
  // test_2d_box = testTreeWithBox<2>();
  // std::cout << "------------- NEXT TEST -------------" << std::endl;
  // test_3d_box = testTreeWithBox<3>();

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  if (lsInternal::all(test_2d_sphere, test_3d_sphere,
                      test_2d_box, test_3d_box))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
};
