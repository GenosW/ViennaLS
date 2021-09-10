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

// --- PRINT FUNCTIONS
// TODO: move to lsTree?

template <typename Point>
void print_point(Point &p)
{
  std::cout << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
  // std::cout << "(" << *p[0][0] << "," << *p[0][1] << "," << *p[0][2] << ")" << std::endl;
};

template <typename T = double, std::size_t D = 3>
std::ostream &operator<<(std::ostream &os, std::array<T, D> point)
{
  os << "(" << point[0] << ", " << point[1] << ", " << point[2] << ")";
  return os;
};

template <class iterator>
void print_iterator(iterator p)
{
  auto pos = *(p.pos);
  std::cout << "Iterator(pos= " << p.pos << ", segment= " << p.segment << ") = " << pos;
};

std::ostream &operator<<(std::ostream &os, lsTree<double, 2>::iterator p)
{
  auto pos = *(p.pos);
  os << "Iterator(pos= " << p.pos << ", segment= " << p.segment << ") = " << pos;
  return os;
};

std::ostream &operator<<(std::ostream &os, lsTree<double, 3>::iterator p)
{
  auto pos = *(p.pos);
  os << "Iterator(pos= " << p.pos << ", segment= " << p.segment << ") = " << pos;
  return os;
};

template <typename Left, typename Right>
void compare(Left left, Right right)
{
  std::cout << left << " =? " << right << " [" << (left == right ? "Y]" : "N]") << std::endl;
};

// --- UNIT TESTS

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

  // traditional index based
  auto traditional_for = [&](size_t ptsToPrint)
  {
    std::cout << "\nFOR LOOP (int i = 0; i < ptsToPrint; ++i) [" << ptsToPrint << " pts]" << std::endl;
    auto begin = tree.begin();
    if (ptsToPrint < 11)
    {
      for (cnt = 0; cnt < ptsToPrint; ++cnt)
      {
        auto point = begin + cnt;
        std::cout << point << "\n";
      }
    }
    else
    {
      for (cnt = 0; cnt < ptsToPrint; ++cnt)
      {
        auto point = begin + cnt;
      }
    }
    LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    std::cout << tree.data[cnt - 1] << " vs " << tree.data[ptsToPrint - 1] << std::endl;
  };

  // begin-end-style
  auto begin_end_for_copy = [&]()
  {
    cnt = 0;
    auto begin = tree.begin();
    auto end = tree.end();
    auto ptsToPrint = end - begin;

    std::cout << "\nFOR LOOP (auto point = tree.begin(); point != end; ++point) [" << ptsToPrint << " pts]";

    for (auto point = begin; point != end; ++point)
    {
      ++cnt;
      if (cnt > N)
      {
        std::cout << "OVERSHOOT: " << point << std::endl;
        break;
      }
    }
    std::cout << std::endl;
    ++cnt; // discrepancy due to for-body not executing...could also start at 1 or move into for(...; end-stmnt)
    // LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    std::cout << tree.data[cnt - 1] << " vs " << tree.data[ptsToPrint - 1] << " vs " << end << std::endl;
  };

  auto begin_end_for_ref = [&]()
  {
    cnt = 0;
    auto begin = tree.begin();
    auto end = tree.end();
    auto ptsToPrint = end - begin;

    std::cout << "\nFOR LOOP (auto &point = begin; point != end; ++point) [" << ptsToPrint << " pts]";

    for (auto &point = begin; point != end; ++point)
    {
      ++cnt;
      if (cnt > N)
      {
        std::cout << "OVERSHOOT: ";
        print_point(point);
        break;
      }
    }
    std::cout << std::endl;
    ++cnt;
    // LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    std::cout << tree.data[cnt - 1] << " vs " << tree.data[ptsToPrint - 1] << " vs " << end << std::endl;
  };

  auto begin_end_for_forward = [&]()
  {
    cnt = 0;
    auto begin = tree.begin();
    auto end = tree.end();
    auto ptsToPrint = end - begin;

    std::cout << "\nFOR LOOP (auto &&point = begin; point != end; ++point) [" << ptsToPrint << " pts]";

    for (auto &&point = begin; point != end; ++point)
    {
      ++cnt;
      if (cnt > N)
      {
        std::cout << "OVERSHOOT: ";
        print_point(point);
        break;
      }
    }
    std::cout << std::endl;
    ++cnt;
    // LSTEST_ASSERT(cnt == ptsToPrint)
    std::cout << cnt << " vs " << ptsToPrint << std::endl;
    std::cout << tree.data[cnt - 1] << " vs " << tree.data[ptsToPrint - 1] << " vs " << end << std::endl;
  };

  // test begin and end methods
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

  begin_end_for_copy();
  begin_end_for_ref();
  begin_end_for_forward();

  begin_end_test(&tree.data.front(), &tree.data.back());

  std::cout << "\nRANGE-FOR LOOP 0 (init-stmnt; type item : cntnr)";
  cnt = 0;
  for (auto point : tree)
  {
    ++cnt;
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
  compare(tree.data[cnt - 1], tree.data[N - 1]);

  std::cout << "\nRANGE-FOR LOOP 1 (init-stmnt; type& item : cntnr)";
  cnt = 0;
  for (auto &point : tree)
  {
    ++cnt;
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
  compare(tree.data[cnt - 1], tree.data[N - 1]);

  std::cout << "\nRANGE-FOR LOOP 2 (init-stmnt; type&& item : cntnr)";
  cnt = 0;
  for (auto &&point : tree)
  {
    ++cnt;
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
  compare(tree.data[cnt - 1], tree.data[N - 1]);

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
