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
  using index_vector = std::vector<size_t>;

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
    bool isLeaf = true;

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

    void leafCheck(size_t colorAssign)
    {
      isLeaf = (left != nullptr || right != nullptr) ? false : true;
      color = colorAssign;
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
  using nodes_vector = std::vector<lsSmartPointer<node>>;
  lsSmartPointer<node> root = nullptr;
  nodes_vector nodes;
  index_vector sortedPoints;

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

    auto nodesBegin = mesh->getNodes().begin();
    auto nodesEnd = mesh->getNodes().end();

    std::vector<lsSmartPointer<index_vector>> orders;

    // 1.Generate absolute orders in each dimension
    // Sort in z/y axis (D=3/D=2)
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    // Sort in x dimension (always)
    size_t dimToSort = 0;
    auto x_idx = lsSmartPointer<index_vector>::New(argsortInDimension(nodesBegin, dimToSort, N));
    orders.push_back(x_idx);

    // Sort in y dimension (if D = 3)
    if constexpr (D > 2)
    {
      dimToSort = 1; // y
      auto y_idx = lsSmartPointer<index_vector>::New(argsortInDimension(nodesBegin, dimToSort, N));
      orders.push_back(y_idx);
    }

    dimToSort = 2; // z
    auto z_idx = lsSmartPointer<index_vector>::New(N);
    std::generate(z_idx->begin(), z_idx->end(), [n = 0]() mutable
                  { return n++; });
    orders.push_back(z_idx);
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    sortedPoints = *z_idx;

    // The absolute, global sorted index-vectors are used to allow for by dimension sorting on node level
    // Example
    // Global order
    // x: 0, 5, 3, 2, 1, 6, 4, 7,
    // y: 5, 1, 0, 4, 6, 3, 7, 2
    // z: 7, 3, 5, 2, 1, 4, 0, 6,
    // kd-Tree by level
    // L1(by x): sortedList = [0, 5,  3, 2] [1, 6,  4, 7]
    // L2(by y): sortedList = [5, 0] [3, 2] [1, 4] [6, 7]
    // L3(by z): sortedList = [5][0] [3][2] [1][4] [7][6]

    size_t index = 0;
    size_t startLastLevel = 0;
    for (size_t level = 0; level < maxDepth + 1; ++level)
    {
      // this levels split dimension
      // 3D:
      // cycle: z -> y -> x -> z -> y -> ...
      // cycle: 2 -> 1 -> 0 -> 2 -> 1 -> ...
      // 2D:
      // cycle: y -> x -> y -> ...
      // cycle: 1 -> 0 -> 1 -> ...
      size_t dim = D - level % (D + 1); // offset
      startLastLevel = buildLevel(sortedPoints, orders, level, dim, N);
    }
    numBins = pow2(depth);
    
    // Array that colors based partitioning
    color = std::vector<T>(N);
    for (size_t i = startLastLevel; i < nodes.size(); ++i)
    {
      auto node = nodes[i];
      for (size_t index = node->start; index < node->stop; ++index)
      {
        color[index] = node->color;
      }
    }
    mesh->insertNextScalarData(color, "color");
  }

  /// Builds the tree level-by-level
  ///
  /// 
  /// TODO: Use Iterators (begin/end) instead of start/stop
  template <class Vector, class VectorPointer>
  size_t buildLevel(Vector &data, VectorPointer &orders, size_t level, size_t dim, size_t N)
  {
    size_t num_nodes = pow2(level);
    size_t nodeSize = (size_t)std::ceil(N / num_nodes); // points per node in the level
    depth = level;                                      // new level --> set tree depth

    // TODO: Just append to nodes vector directly
    // Just here, cause unsure if needed.
    nodes_vector levelNodes(num_nodes);
    auto start = 0;
    auto stop = start + nodeSize;
    auto order_begin = orders[dim]->begin();
    auto order_end = orders[dim]->end();
    auto data_begin = data.begin();
    auto data_end = data.end();
    for (auto &thisNode : levelNodes)
    {
      thisNode = lsSmartPointer<node>::New();
      thisNode->setRange(start, stop);
      sortByIdxRange(data_begin + start, data_begin + stop, order_begin, order_end);
      thisNode->level = level;
      thisNode->dimSplit = dim;

      // TODO: Connect tree
    }

    bool isFinal = ((nodeSize < maxPointsPerBin) || (level > maxDepth)) ? false : true;
    if (isFinal) {
      int col = -int(levelNodes.size() / 2);
      for (auto &thisNode : levelNodes)
      {
        thisNode->color = col;
        ++col;
      }
    }
    size_t thisLevelStart = nodes.size();
    nodes.insert(nodes.end(), levelNodes.begin(), levelNodes.end());
    return thisLevelStart;
  }

  /// builds a node of the tree and returns a pointer
  ///
  /// checks if it should be a leaf node
  /// TODO: Use Iterators (begin/end) instead of start/stop
  template <class VectorIt>
  lsSmartPointer<node> buildByDepth(VectorIt begin, size_t start, size_t stop, size_t level, size_t &index)
  {
    size_t size = stop - start;
    if ((size < maxPointsPerBin) || (level > maxDepth))
      return nullptr;
    size_t halfSize = size / 2;
    depth = (depth < level) ? level : depth;
    //size_t thisIndex = index++;

    // thisNodes split dimension
    size_t dim = D - index % (D + 1); // offset
    ++index;

    auto thisNode = lsSmartPointer<node>::New();
    nodes.push_back(thisNode);
    thisNode->setRange(start, stop);
    thisNode->level = level;
    thisNode->dimSplit = dim;

    thisNode->left = build(begin, start, start + halfSize, level + 1, index);
    thisNode->right = build(begin, start + halfSize, stop, level + 1, index);
    //thisNode->isLeaf = (thisNode->left != nullptr || thisNode->right != nullptr ) ? false : true;
    thisNode->leafCheck(getNextLeafColor());

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

  /// Sorts a vector (range) by the Reihenfolge given by second vector
  /// Elements of the first, that are not present in the second,
  /// are pushed to the end.
  ///
  /// Example
  /// x: 0, 5, 3, 2, 1, 6, 4, 7,
  /// y: 5, 1, 0, 4, 6, 3, 7, 2
  /// z: 7, 3, 5, 2, 1, 4, 0, 6,
  /// L1(by x): sortedList = [0, 5,  3, 2] [1, 6,  4, 7]
  /// L2(by y): sortedList = [5, 0] [3, 2] [1, 4] [6, 7]
  /// L3(by z): sortedList = [5][0] [3][2] [1][4] [7][6]
  template <class It, class It2>
  void sortByIdxRange(It begin, It end, It2 beginSortBy, It2 endSortBy)
  {
    std::sort(begin, end,
              [beginSortBy, endSortBy](auto left, auto right)
              {
                // sort in ascending order
                for (auto it = beginSortBy; it < endSortBy; ++it)
                {
                  if (*it == left)
                    return true;
                  if (*it == right)
                    return false;
                }
                return false;
              });
  }

  /// Sorts a vector of 3D points (point cloud) by generating a sorted
  /// vector of indices (of length size) and move-returns it.
  /// Does NOT actually sort the vector (in place)!
  template <class VectorIt>
  index_vector argsortInDimension(VectorIt toSortBegin, size_t dim, size_t size)
  {
    std::vector<size_t> result(size);
    std::generate(result.begin(), result.end(), [n = 0]() mutable
                  { return n++; });
    // Taken and modified from: https://stackoverflow.com/questions/2312737/whats-a-good-way-of-temporarily-sorting-a-vector
    std::sort(result.begin(), result.end(),
              [toSortBegin, dim](auto left, auto right) {                // left, right : size_t -> indices of vector toSort
                return toSortBegin[left][dim] < toSortBegin[right][dim]; // sort in ascending order
              });
    return std::move(result);
  }

  /// OLD sortInDim
  /// TODO: Remove old sort
  template <class VectorIt, class VectorItRes>
  index_vector partitionInDimensionOLD(VectorIt toSortBegin, VectorItRes resultBegin, VectorItRes resultEnd, size_t dim)
  {
    // Taken and modified from: https://stackoverflow.com/questions/2312737/whats-a-good-way-of-temporarily-sorting-a-vector
    std::sort(resultBegin, resultEnd,
              [toSortBegin, dim](auto left, auto right) {                // left, right : size_t -> indices
                return toSortBegin[left][dim] < toSortBegin[right][dim]; // sort in ascending order
              });
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
  * ---- UTILITY METHODS ---- *
  */
private:
  size_t getNextLeafColor()
  {
    return 0;
  }

  /*
  * ---- DEBUG METHODS ---- *
  */
public:
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
