//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdint>
#include <lsToDiskMesh.hpp>
#include <lsSmartPointer.hpp>

constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1)
{
  return exp < 1 ? result : ipow(base * base, exp / 2, (exp % 2) ? result * base : result);
}

// all depth values are known at compile time
constexpr int64_t pow2(int exp)
{
  return 1 << exp;
}

inline int ipower(int N, int exp) { return (exp > 1) ? N * ipower(N, exp - 1) : N; };

template <class Type>
struct SqrDiff
{
  Type operator()(Type x, Type y) const
  {
    return (x - y) * (x - y);
  }
};

template <typename T, class Iter_T, class Iter2_T>
T pointDistance(Iter_T first, Iter_T last, Iter2_T first2)
{
  T ret = inner_product(first, last, first2, 0.0,
                        std::plus<T>(), SqrDiff<T>());
  return ret > 0.0 ? sqrt(ret) : 0.0;
}

/// Tree structure that seperates a domain
///
/// TODO: clean up parameters, member functions and attributes
///
/// Maximum depth is 4?
template <class T, int D>
class lsTree
{
public:
  using point_type = std::array<T, 3>;
  // typedef typename std::array<T, 3> point_type2;

private:
  // The mesh we're building a tree for
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  // lsSmartPointer<lsMesh<T>> newMesh = nullptr;
  std::vector<size_t> x;
  std::vector<size_t> y;
  std::vector<size_t> z;
  std::vector<T> color;

  // PARAMETERS
  std::string tree_type = "kd-Tree";
  int depth = 0;
  int numBins = 1; // numBins = 2^depth
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  int maxDepth = 4;
  int maxNumBins = pow2(maxDepth);
  int maxPointsPerBin = 10;

  struct node
  {
    // std::vector<point_type> points;
    size_t start = 0;
    size_t stop = 0;
    size_t level = 0;
    size_t size = 0;
    size_t color = 0;
    size_t dimSplit = 0;
    lsSmartPointer<node> left = nullptr;
    lsSmartPointer<node> right = nullptr;

    point_type center; // unused for now
    point_type minimumExtent;
    point_type maximumExtent;

    // node(const point_type& ct, const point_type& tl, const point_type& br) : center(ct), extent({tl, br}) left(nullptr), right(nullptr) {
    // }
    node()
    {
    }

    node(const point_type &ct) : center(ct)
    {
    }

    void setRange(size_t newStart, size_t newStop)
    {
      start = newStart;
      stop = newStop;
      size = stop - start;
    }

    void setCenter(point_type ct_point)
    {
      center = ct_point;
    }

    point_type &getCenter() const
    {
      return center;
    }

    bool isWithin(const point_type &pt)
    {
      for (size_t i = 0; i < 3; ++i)
        if (pt[i] < minimumExtent || pt[i] > maximumExtent)
          return false;
      return true;
    }

    T distance(const point_type &pt)
    {
      return pointDistance<T>(center.begin(), center.end(), pt.begin());
    }

    /// Unused (from discarded attempt)
    // void split()
    // {
    //   left = lsSmartPointer<node>::New();
    //   size_t med = (stop - start) / 2;
    //   left->setRange(start, med);
    //   right = lsSmartPointer<node>::New();
    //   right->setRange(med, stop);
    // }
  };
  lsSmartPointer<node> root = nullptr;
  std::vector<lsSmartPointer<node>> nodes;

public:
  lsTree() {}

  lsTree(lsSmartPointer<lsMesh<T>> passedMesh)
      : mesh(passedMesh) {}

  void setMesh(lsSmartPointer<lsMesh<T>> passedMesh)
  {
    mesh = passedMesh;
  }

  //~lsTree() {}

  /// entry point
  void apply()
  {
    if (mesh == nullptr)
    {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsTree.")
          .print();
      return;
    }
    nodes.reserve(maxNumBins);
    // partition in x by median
    size_t N = mesh->getNodes().size();
    size_t medianPos = size_t(N / 2);

    std::cout << "N: " << N << std::endl;
    std::cout << "medianPos: " << medianPos << std::endl;

    auto begin = mesh->getNodes().begin();
    auto end = mesh->getNodes().end();

    // Partition in z/y axis
    size_t dimToPart = 2; // z
    // Array that colors based on z (D=3)/y (D=2) partitioning (just for output)
    color = std::vector<T>(N, -1);
    std::for_each(color.begin() + medianPos, color.end(), [](auto &item) { item = 1; });

    // Partition other dimensions
    if (D > 2)
    {
      dimToPart = 1; // y
      y = std::vector<size_t>(N, 0);
      std::generate(y.begin(), y.end(), [n = 0]() mutable { return n++; });
      partitionInDimension(begin, y.begin(), y.end(), dimToPart));
    }
    dimToPart = 0;
    x = std::vector<size_t>(N);
    n = 0;
    std::generate(x.begin(), x.end(), [n = 0]() mutable { return n++; });
    partitionInDimension(begin, x.begin(), x.end(), dimToPart));

    size_t index = 0;
    root = build(begin, 0, N, 0, index);
    numBins = pow2(depth);
    mesh->insertNextScalarData(color, "color");
  }

  /// builds a node of the tree and returns a pointer
  ///
  /// checks if it should be a leaf node
  /// TODO: add partitioning of vector/mesh-nodes
  /// TODO: Use Iterators (begin/end) instead of start/stop
  template <class VectorIt>
  lsSmartPointer<node> build(VectorIt begin, size_t start, size_t stop, size_t level, size_t &index)
  {
    size_t size = stop - start;
    if ((size < maxPointsPerBin) || (level > maxDepth))
      return nullptr;
    size_t halfSize = size / 2;
    // TODO: add partitioning of vector/mesh-nodes
    //        or use partitioning here
    depth = (depth < level) ? level : depth;
    //size_t thisIndex = index++;

    // thisNodes split dimension
    size_t dim = D - index++ % (D + 1); // offset

    auto thisNode = lsSmartPointer<node>::New();
    nodes.push_back(thisNode);
    thisNode->setRange(start, stop);
    thisNode->level = level;
    thisNode->dim = dim;

    thisNode->left = build(begin, start, start + halfSize, level + 1, index);
    thisNode->right = build(begin, start + halfSize, stop, level + 1, index);

    return thisNode;
  }

  /// [WIP] query the tree based on a single point
  ///
  /// checks if it should be a leaf node
  /// TODO: TEST
  const lsSmartPointer<node> nearest(const point_type &pt)
  {
    if (root == nullptr)
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning pt.")
          .print();
      return nullptr;
    }

    lsSmartPointer<node> thisNode = root;
    while (thisNode->left != nullptr)
      thisNode = (thisNode->left->isWithin(pt)) ? thisNode->left : thisNode->right;
    return thisNode;
  }

  template <class VectorIt, class VectorItRes>
  size_t partitionInDimension(VectorIt toSort, VectorItRes resultBegin, VectorItRes resultEnd, size_t dim)
  {
    // Taken and modified from: https://stackoverflow.com/questions/2312737/whats-a-good-way-of-temporarily-sorting-a-vector
    std::sort(resultBegin, resultEnd,
              [toSort](int left, int right) {
                return toSort[left][dim] < toSort[right][dim]; // sort in ascending order
              });

    // previous
    // std::partition(begin, end, [dim, median](const auto &pos) { return pos[dim] < median; }) return std::distance(begin, end) / 2
  }

  // bool isWithinExtent(point_type &pt)
  // {
  //   for (size_t i = 0; i < pt.size(); ++i)
  //   {
  //     if (pt[i] < mesh->minimumExtent[i] || pt[i] > mesh->maximumExtent[i])
  //       return false;
  //   }
  //   return true;
  // }

  /*
  * ---- DEBUG FUNCTIONS ---- *
  */
  const std::string getTreeType()
  {
    return tree_type;
  }

  /// prints parameters
  void printInfo()
  {
    std::cout << getTreeType() << std::endl;
    std::cout << "depth: " << depth << std::endl;
    std::cout << "numBins: " << numBins << std::endl;
    std::cout << "maxDepth: " << maxDepth << std::endl;
    std::cout << "maxNumBins: " << maxNumBins << std::endl;
    std::cout << "maxPointsPerBin: " << maxPointsPerBin << std::endl;
  };

  /// prints tree by simply going through nodes-vector
  void printTree()
  {
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      std::string leaf = (nodes[i]->left == nullptr && nodes[i]->right == nullptr) ? "#" : "";
      std::cout << "node " << std::setw(3) << i << "(L" << nodes[i]->level << "): [" << nodes[i]->start << ", " << nodes[i]->stop << ")" << leaf << std::endl;
    }
  };

  /// prints tree one level per line at a time
  void printTreeByLevel()
  {
    for (size_t level = 0; level <= depth; ++level)
    {
      std::cout << "L" << level << ": ";
      for (size_t i = 0; i < nodes.size(); ++i)
      {
        if (nodes[i]->level != level)
          continue;
        std::cout << std::setw(3) << i << "[" << nodes[i]->start << ", " << nodes[i]->stop << ") --- ";
      }
      std::cout << std::endl;
    }
  };
};

#endif // LS_TREE
