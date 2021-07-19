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
#include <lsTestAsserts.hpp>
//#define LS_TREE
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>

#define PRINT(Var) std::cout << #Var << ": " << Var

enum struct lsTestStatus : unsigned
{
  SUCCESS = 0,
  FAILED = 1,
  UNCHECKED = 2

};

struct lsTest
{
  lsTestStatus status = lsTestStatus::UNCHECKED;
  std::string name = "unnamed";

  lsTest(const lsTestStatus stat)
  {
    status = stat;
  }
  lsTest(lsTestStatus stat, std::string nam) : status(stat), name(nam) {}

  template <typename T>
  bool run(T (&test)(void))
  {
    this->status = lsTestStatus{(uint)test()};
    return this->status == lsTestStatus::SUCCESS;
  }

  void check() const
  {
    std::cout << name;
    if (status == lsTestStatus::FAILED)
    {
      std::cout << ": FAILED";
    }
    else if (status == lsTestStatus::SUCCESS)
    {
      std::cout << ": SUCCESS";
    }
    else if (status == lsTestStatus::UNCHECKED)
    {
      std::cout << ": UNCHECKED";
    }
    else
    {
      std::cout << ": unknown status";
    }
    std::cout << std::endl;
  };

  lsTestStatus operator()() const
  {
    return this->status;
  }

  bool wasSuccess() const
  {
    return this->status == lsTestStatus::SUCCESS;
  }
};

#define INIT_LSTEST(Var) lsTest(lsTestStatus::UNCHECKED, #Var)
#define MAKE_LSTEST(Var) lsTest Var(lsTestStatus::UNCHECKED, #Var)

template <class Test, class... Tests>
void check(Test const &t1, Tests const &...rest)
{
  t1.check();
  if constexpr (sizeof...(rest) > 0)
  {
    check(rest...);
  }
}

template <typename... Args>
bool all(Args... args)
{
  return (... && args);
}

template <typename T, typename... Args>
bool all_equal(T target, Args... args)
{
  return ((target == args) && ...);
}

void setup_omp(uint numThreads)
{
  int maxThreads = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << maxThreads << ") threads for this test." << std::endl;
};

template <int Dim>
lsTestStatus testTree(void)
{
  constexpr int D = Dim;
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

  // Get reference geometry parameters
  size_t N = levelSet->getDomain().getNumberOfPoints();
  std::cout << "Initial: " << std::endl;
  std::cout << "Number of points: " << N << std::endl;

  // --- Volumetric mesh
  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, prefix + "_Mesh.vtk").apply();

  // Compute lsTree of volumetric mesh
  auto tree = lsTree<double, D>(mesh);
  tree.apply();
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh.vtu").apply();

  tree.printBT();
  // EVALUATE
  auto numLevels = tree.getDepth();
  std::vector<size_t> ptsInTree(numLevels + 1, 0);
  for (auto &node : tree.getTreeNodes())
  {
    // std::cout << "{color: " << node->color
    //           << " / level: " << node->level
    //           << " / isLeaf: " << node->isLeaf
    //           << "} -> size:" << node->size() << std::endl;
    ptsInTree[node->level] += node->size();
  }
  LSTEST_ASSERT(
      std::all_of(ptsInTree.cbegin(), ptsInTree.cend(),
                  [N](auto item)
                  { return item == N; }));

  // --- DiskMesh
  // std::cout << "DiskMesh: " << std::endl;
  // lsToDiskMesh<double, D>(levelSet, mesh).apply();
  // lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_DiskMesh.vtu").apply();

  // auto tree2 = lsTree<double, D>(mesh);
  // tree2.apply();
  // lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeDiskMesh.vtu").apply();

  // numLevels = tree2.getDepth();
  // ptsInTree = std::vector<size_t>(numLevels + 1, 0);
  // for (auto &node : tree.getTreeNodes())
  // {
  //   // std::cout << "{color: " << node->color
  //   //           << " / level: " << node->level
  //   //           << " / isLeaf: " << node->isLeaf
  //   //           << "} -> size:" << node->size() << std::endl;
  //   ptsInTree[node->level] += node->size();
  // }
  // LSTEST_ASSERT(
  //     std::all_of(ptsInTree.cbegin(), ptsInTree.cend(),
  //                 [N](auto item)
  //                 { return item == N; }));

  // PASS!
  return lsTestStatus::SUCCESS;
};

int main(int argc, char **argv)
{
  lsTest test_2d = INIT_LSTEST(test_2d);
  MAKE_LSTEST(test_3d);

  test_2d.run(testTree<2>);
  // test_3d.run(testTree<3>);

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  check(test_2d, test_3d);

  if (all(test_2d.wasSuccess(), test_3d.wasSuccess()))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
}
