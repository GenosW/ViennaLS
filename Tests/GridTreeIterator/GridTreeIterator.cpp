// STL
#include <iostream>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsDimension.hpp>
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>
#include <lsTest.hpp>
#include <lsTestAsserts.hpp>

template <typename T = double, int D = 3>
std::ostream &operator<<(std::ostream &os, std::array<T, D> point)
{
  os << "(" << point[0] << ", " << point[1] << ", " << point[2] << ")";
  return os;
};

// template <typename T = double, int D = 3>
// std::ostream &operator<<(std::ostream &os, typename lsTree<T, D>::template treeIterator<std::array<T, D>> point)
// {
//   os << "(" << point[0] << ", " << point[1] << ", " << point[2] << ")";
//   return os;
// };

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

template <int D>
uint testTreeIterator(void)
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
  checkTree(tree, N);
  lsVTKWriter(mesh, lsFileFormatEnum::VTU, prefix + "_TreeMesh").apply();

  // std::array<double, 3> point{bounds[0] + 1, bounds[2] + 1, bounds[4]};
  // auto treeIt = tree2.getNearestNeighbors(point);
  // std::cout << *treeIt;

  // Test treeIterator
  int cnt;
  auto print_iterator = [](auto p)
  {
    std::cout << "Iterator(pos= " << p.pos << ", segment= " << p.segment << ")"
              << std::endl;
    auto pos = *(p.pos);
    std::cout << "(" << pos[0] << ")" << std::endl;
  };

  auto print_point_pointer = [](auto p)
  {
    std::cout << "(" << (*p)[0] << ", " << p[1] << ", " << p[2] << ")" << std::endl;
  };

  auto print_point = [](auto p)
  {
    std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
    // std::cout << "(" << *p[0][0] << "," << *p[0][1] << "," << *p[0][2] << ")" << std::endl;
  };

  auto compare = [](auto left, auto right)
  {
    std::cout << left << " =? " << right << " [" << (left == right ? "Y]" : "N]") << std::endl;
  };

  auto traditional_for = [&](size_t ptsToPrint)
  {
    cnt = 0;
    std::cout << "\nFOR LOOP (traditional) [" << ptsToPrint << " pts]" << std::endl;
    // for (auto point : tree)
    for (auto point = tree.begin(); cnt < ptsToPrint; ++point)
    // for (auto point = tree.data.begin(); cnt < ptsToPrint; ++point)
    {
      ++cnt;
      // print_point_pointer(point);
    }
    // std::cout << std::endl;
    LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    print_point(tree.data[cnt - 1]);
    std::cout << " vs ";
    print_point(tree.data[ptsToPrint - 1]);
    std::cout << std::endl;
  };

  auto begin_end_for = [&](size_t ptsToPrint)
  {
    cnt = 0;
    std::cout << "\nFOR LOOP (begin-end) [" << ptsToPrint << " pts]";
    // for (auto point : tree)
    auto end = tree.end();
    for (auto &&point = tree.begin(); point != end; ++point)
    // for (auto point = tree.data.begin(); cnt < ptsToPrint; ++point)
    {
      ++cnt;
      if (cnt > N)
      {
        std::cout << "OVERSHOOT: ";
        print_point_pointer(point);
        break;
      }
    }
    std::cout << std::endl;
    ++cnt;
    // LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    // compare()
    print_point(tree.data[cnt - 1]);
    std::cout << " vs ";
    print_point(tree.data[ptsToPrint - 1]);
    std::cout << " vs ";
    print_point(end);
    std::cout << std::endl;
  };

  auto begin_end_test = [&](auto beginShould, auto endShould)
  {
    std::cout << "\nBEGIN" << std::endl;
    auto begin = tree.begin();
    auto im_lazy = [&](auto left, auto right)
    {
      print_point(left);
      std::cout << " =? ";
      print_point(right);
      std::cout << " [" << (*left == right ? "Y]" : " N]") << std::endl;
    };
    im_lazy(begin, tree.data[0]);
    compare(begin.pos, beginShould);

    std::cout << "END" << std::endl;
    auto end = tree.end();
    im_lazy(end, tree.data.back());
    compare(end.pos, endShould);
  };

  traditional_for(10);
  traditional_for(N);
  begin_end_for(N);
  begin_end_test(&tree.data.front(), &tree.data.back());

  std::cout << "\nRANGE-FOR LOOP (init-stmnt; type item : cntnr)";
  cnt = 0;
  for (auto &point : tree)
  {
    ++cnt;
    // std::cout << ++cnt << ": (" << point[0] << ", "
    //           << point[1] << ", "
    //           << point[2] << ", "
    //           << ")\n";
    if (cnt == N)
    {
      std::cout << "OVERSHOOT: " << &point << " aka ("
                << point[0] << ", "
                << point[1] << ", "
                << point[2] << ")"
                << std::endl;
      break;
    }
  }
  ++cnt;
  std::cout << std::endl;
  // LSTEST_ASSERT(cnt == N)
  std::cout << cnt << " vs " << N << std::endl;
  print_point(tree.data[cnt - 1]);
  std::cout << " vs ";
  print_point(tree.data[N - 1]);
  std::cout << std::endl;

  std::cout << "\nRANGE-FOR LOOP on std::vector" << std::endl;
  std::vector<int> testitest(10, -1);
  int cnt2 = 0;
  for (auto val : testitest)
    ++cnt2;
  std::cout << cnt2 << std::endl;

  // PASS!
  return 1;
};

int main(int argc, char **argv)
{
  uint test_2d_box = 0;
  uint test_3d_box = 1;

  std::cout << "################ BOX ################" << std::endl;
  test_2d_box = testTreeIterator<2>();
  std::cout << "------------- NEXT TEST -------------" << std::endl;
  // test_3d_box = testTreeIterator<3>();

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  if (lsInternal::all(test_2d_box, test_3d_box))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return EXIT_SUCCESS;
};
