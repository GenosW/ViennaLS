//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <iterator>
#include <vector>
#include <cstdint>
#include <lsToDiskMesh.hpp>
#include <lsSmartPointer.hpp>

#pragma region helpers
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
#pragma endregion
/// Tree structure that seperates a domain
///
/// TODO: clean up parameters, member functions and attributes
///
/// Maximum depth is 4?
template <class T, int D>
class lsTree
{
public:
  using value_type = T;
  using size_type = std::size_t;

  using point_type = std::array<T, 3>;
  using index_vector = std::vector<size_t>;
  using data_vector = std::vector<point_type>; // TODO: Template this

private:
  // The mesh we're building a tree for
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  //lsSmartPointer<std::vector<T>> dataToSort = nullptr;
  // lsSmartPointer<lsMesh<T>> newMesh = nullptr;
  index_vector x;
  index_vector y;
  index_vector z;
  std::vector<T> color;

  // PARAMETERS
  std::string tree_type = "kd-Tree";
  int depth = 0;
  int numBins = 1; // numBins = 2^depth
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  int maxDepth = 4;
  int maxNumBins = pow2(maxDepth);
  int maxPointsPerBin = 4;

#pragma region treeNode
  struct treeNode
  {
    bool isLeaf = false;
    size_t start = 0;
    size_t stop = 0;
    size_t level = 0;
    //size_t size = 0;
    int color = 0;
    size_t dimSplit = 0;

    lsSmartPointer<treeNode> left = nullptr;
    lsSmartPointer<treeNode> right = nullptr;

    T median;

    // TODO: Evaluate these
    point_type center; // unused for now

    // treeNode(const point_type& ct, const point_type& tl, const point_type& br) : center(ct), extent({tl, br}) left(nullptr), right(nullptr) {
    // }
    treeNode()
    {
    }

    treeNode(size_t lev, size_t dim, size_t startRange, size_t stopRange, T med) : level(lev), dimSplit(dim), start(startRange), stop(stopRange), median(med)
    {
    }

    void setRange(size_t newStart, size_t newStop)
    {
      start = newStart;
      stop = newStop;
      //size = stop - start;
    }

    size_t size()
    {
      return stop - start;
    }

    bool belowMedian(const point_type &pt)
    {
      // TODO: implement
      return pt[dimSplit] <= median;
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

    void setLeaf(bool value)
    {
      isLeaf = value;
    }

    /// Unused (from discarded attempt)
    // void split()
    // {
    //   left = lsSmartPointer<treeNode>::New();
    //   size_t med = (stop - start) / 2;
    //   left->setRange(start, med);
    //   right = lsSmartPointer<treeNode>::New();
    //   right->setRange(med, stop);
    // }
  };
#pragma endregion
  using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
  nodes_vector treeNodes;
  index_vector sortedPoints;
  data_vector

#pragma region treeIterator
      template <class T_base = T>
      struct treeIterator
  {
  public:
    using iterator_category = std::forward_iterator_tag; //std::bidirectional_iterator_tag;
    using value_type = T_base;
    using difference_type = size_type;
    using pointer = value_type *;
    using reference = value_type &;

    std::vector<std::array<size_t, 2>> ranges;
    // Idea: Utilize [start, stop) of each leaf node
    // ranges = [(startBin1, stopBin1), (startBin2, stopBin2),
    //           (startBin3, stopBin3), (startBin4, stopBin4), ... ]
    size_t pos;

  public:
    treeIterator() = default;
    treeIterator(treeNode *init) : pos(init->start){}; // probably superfluous...
    // treeIterator(const MyConstIterator& other) = delete; // no need to delete...doesnt exist anyway

    reference operator*()
    {
      return this->pos->value;
    }

    pointer operator->()
    {
      return &(this->pos->value);
    }

    bool operator==(const treeIterator &other)
    {
      return this->pos == other.pos;
    }

    bool operator!=(const treeIterator &other)
    {
      return this->pos != other.pos;
    }

    treeIterator &operator++()
    {
      this->pos = this->pos->next;
      return *this;
    }

    treeIterator &operator--()
    {
      this->pos = this->pos->previous;
      return *this;
    }
  };

  template <class T_it>
  struct MyConstIterator : treeIterator<const T_it>
  {
    // This is all not needed and actually problematic (shadowing...)
    // public:
    // using iterator_category = std::bidirectional_iterator_tag;
    // using difference_type   = size_type; // std::ptrdiff_t;
    // using value_type        = const T_it;
    // using pointer           = value_type*;  // or also value_type*
    // using reference         = value_type&;  // or also value_type&

    // Node* pos = nullptr;

  public:
    // MyConstIterator() = default; // default constructor enabled
    using treeIterator<const T_it>::treeIterator;

    /** Construct from non-const iterator
  */
    MyConstIterator(const treeIterator<T_it> &nonconst_other)
    {
      this->pos = nonconst_other.pos;
    };
  };

  using iterator = treeIterator<T>;
  using const_iterator = MyConstIterator<T>;

  // methods using iteratos
  iterator begin()
  {
    return iterator{this->dataHead};
  };
  iterator end()
  {
    return iterator{this->pastTheEnd};
  }
  const_iterator begin() const
  {
    return const_iterator{this->dataHead};
  };
  const_iterator end() const
  {
    return const_iterator{this->pastTheEnd};
  };

  const_iterator cbegin() const
  {
    return const_iterator{this->dataHead};
  };
  const_iterator cend() const
  {
    return const_iterator{this->pastTheEnd};
  };
#pragma endregion

private:
public:
  lsTree()
  {
  }

  lsTree(lsSmartPointer<lsMesh<T>> passedMesh)
      : mesh(passedMesh) {}

  void setMesh(lsSmartPointer<lsMesh<T>> passedMesh)
  {
    mesh = passedMesh;
  }

  //~lsTree() {}

  /// entry point
  void apply(bool byLevel = true)
  {
    if (mesh == nullptr)
    {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsTree.")
          .print();
      return;
    }
    //treeNodes.reserve(maxNumBins);
    // partition in x by median
    size_t N = mesh->getNodes().size();
    size_t medianPos = size_t(N / 2);

    std::cout << "N: " << N << std::endl;
    std::cout << "medianPos: " << medianPos << std::endl;

    auto dataNodes = mesh->getNodes();

    std::vector<lsSmartPointer<index_vector>> orders;

    // 1.Generate absolute orders in each dimension
    // Sort in z/y axis (D=3/D=2)
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    // Sort in x dimension (always)
    size_t dimToSort = 0;
    auto x_idx = lsSmartPointer<index_vector>::New(argsortInDimension(dataNodes.begin(), dimToSort, N));
    orders.push_back(x_idx);

    // Sort in y dimension (if D = 3)
    if constexpr (D > 2)
    {
      dimToSort++; // dimToSort = 1 <==> y
      auto y_idx = lsSmartPointer<index_vector>::New(largsortInDimension(dataNodes.begin(), dimToSort, N));
      orders.push_back(y_idx);
    }

    dimToSort++; // dimToSort = 2 <==>  z
    auto z_idx = lsSmartPointer<index_vector>::New(N);
    std::generate(z_idx->begin(), z_idx->end(), [n = 0]() mutable
                  { return n++; });
    orders.push_back(z_idx);
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    sortedPoints = index_vector(N);
    std::generate(sortedPoints.begin(), sortedPoints.end(), [n = 0]() mutable
                  { return n++; });

    // The absolute, global sorted index-vectors are used to allow for by dimension sorting on treeNode level
    // Example
    // Global order
    // x: 0, 5, 3, 2, 1, 6, 4, 7,
    // y: 5, 1, 0, 4, 6, 3, 7, 2
    // z: 7, 3, 5, 2, 1, 4, 0, 6,
    // kd-Tree by level
    // L1(by x): sortedList = [0, 5,  3, 2] [1, 6,  4, 7]
    // L2(by y): sortedList = [5, 0] [3, 2] [1, 4] [6, 7]
    // L3(by z): sortedList = [5][0] [3][2] [1][4] [7][6]

    size_t startLastLevel = 0;
    if (byLevel)
    {
      //size_t index = 0;
      size_t binSize = N;
      auto root = lsSmartPointer<treeNode>::New(0, dimToSort, 0, N, dataNodes[sortedPoints[N / 2]][dimToSort]);
      treeNodes.push_back(root);
      for (size_t level = 1; (level < maxDepth + 1) && (binSize > maxPointsPerBin); ++level)
      {
        // this levels split dimension
        // 3D:
        // cycle: z -> y -> x -> z -> y -> ...
        // cycle: 2 -> 1 -> 0 -> 2 -> 1 -> ...
        // 2D:
        // cycle: y -> x -> y -> ...
        // cycle: 1 -> 0 -> 1 -> ...
        std::cout << "D: " << D << std::endl;
        size_t dim = D - 1 - level % (D); // offset
        startLastLevel = buildLevel(sortedPoints, orders, dataNodes, level, dim, N);
        binSize = treeNodes.back()->size();
      }
      numBins = pow2(depth);
    }
    if (!byLevel)
    {
      auto root = buildByDepth(sortedPoints, orders, dataNodes, 0, N, 0);
    }

    if (byLevel)
    {
      size_t idx = 0;
      for (auto &treeNode : treeNodes)
      {
        treeNode->color = idx;
        ++idx;
      }
      // Array that colors based partitioning
      color = std::vector<T>(N);
      for (size_t i = startLastLevel; i < treeNodes.size(); ++i)
      {
        auto treeNode = treeNodes[i];
        for (size_t index = treeNode->start; index < treeNode->stop; ++index)
        {
          //color[index] = treeNode->color;
          color[index] = treeNode->color;
        }
      }
      mesh->insertNextScalarData(color, "color");
    }
  }

  /// Builds the tree level-by-level
  ///
  ///
  /// TODO: Use Iterators (begin/end) instead of start/stop
  template <class Vector, class VectorOfPointers, class VectorData>
  size_t buildLevel(Vector &data, VectorOfPointers &orders, VectorData &originalData, size_t level, size_t dim, size_t N)
  {
    size_t num_nodes = pow2(level);
    // TODO:
    //size_t nodeSizeR = (size_t)(N / num_nodes); // points per treeNode in the level
    depth = level; // new level --> set tree depth
    std::cout << "Level: " << level << " [dim=" << dim << "]" << std::endl;
    std::cout << "num_nodes: " << num_nodes << std::endl;
    std::cout << "nodeSize: ca. " << (N / num_nodes) << std::endl;
    std::cout << "treeNodes.size(): " << treeNodes.size() << std::endl;
    std::cout << "parent nodes: " << treeNodes.size() - (pow2(level - 1) - 1) << std::endl;

    // TODO: Just append to treeNodes vector directly
    // Just here, cause unsure if needed.
    nodes_vector levelNodes;

    auto order_begin = orders[dim]->begin();
    auto order_end = orders[dim]->end();
    auto data_begin = data.begin();
    auto data_end = data.end();

    size_t num_parents = treeNodes.size();

    // TODO: rework to work starting from the previous levels treeNodes as roots and connect that way
    // still save in vector levelNodes for now and then directly append to treeNodes-vector
    // for (auto rootsIt = treeNodes.begin() + pow2(level - 1) - 1; rootsIt < treeNodes.end(); ++rootsIt)
    for (auto idx = pow2(level - 1) - 1; idx < num_parents; ++idx)
    {
      std::cout << "Root " << idx << " - ";
      // root is to be split, so is not a leaf anymore
      // auto root = *rootsIt;
      lsSmartPointer<treeNode> root = treeNodes[idx];
      root->isLeaf = false;

      // Divide into two child nodes
      // Calculate start-, stop-indices and ranges
      size_t start = root->start;
      size_t stop = root->stop;
      size_t range = stop - start;
      size_t leftRange = range / 2;
      size_t rightRange = range - leftRange;
      // Sort Range of the root
      sortByIdxRange(data_begin + start, data_begin + stop, order_begin, order_end);

      // Make 2 new nodes
      // treeNode(level, dim, startIndex, stopIndex, median)
      // median = lookup Point from original data (only )
      T median = originalData[sortedPoints[start + leftRange / 2]][dim];
      lsSmartPointer<treeNode> left = lsSmartPointer<treeNode>::New(level, dim, start, start + leftRange, median);
      //left->median = originalData[sortedPoints[start + leftRange / 2]][dim];
      root->left = left;

      median = originalData[sortedPoints[start + leftRange + rightRange / 2]][dim];
      lsSmartPointer<treeNode> right = lsSmartPointer<treeNode>::New(level, dim, start + leftRange, stop, median);
      //right->median = originalData[sortedPoints[start + leftRange + rightRange / 2]][dim];
      root->right = right;

      // Insert pointers to new nodes into vector
      levelNodes.push_back(left);
      levelNodes.push_back(right);
    }

    //bool isFinal = ((nodeSize < maxPointsPerBin) || (level > maxDepth)) ? true : false;
    // bool isFinal = (level >= maxDepth) ? true : false;
    // if (isFinal)
    // {
    //   int col = -int(levelNodes.size() / 2);
    //   for (auto &thisNode : levelNodes)
    //   {
    //     thisNode->isLeaf = true;
    //     thisNode->color = col;
    //     ++col;
    //   }
    // }
    size_t thisLevelStart = treeNodes.size();
    treeNodes.insert(treeNodes.end(), levelNodes.begin(), levelNodes.end());
    std::cout << "\n----------------- Level " << level << " DONE -----------------" << std::endl;
    return thisLevelStart;
  }

  /// builds a treeNode of the tree and returns a pointer
  ///
  /// checks if it should be a leaf treeNode
  /// TODO: Messy template
  template <class Vector, class VectorOfPointers, class VectorData>
  lsSmartPointer<treeNode> buildByDepth(Vector &data, VectorOfPointers &orders, VectorData &originalData, size_t start, size_t stop, size_t level)
  {
    size_t range = stop - start;
    if ((range < maxPointsPerBin) || (level > maxDepth))
      return nullptr;
    size_t leftRange = range / 2;
    depth = (depth < level) ? level : depth;
    //size_t thisIndex = index++;

    // thisNodes split dimension
    size_t dim = D - level % (D + 1); // offset
    //++index;

    T median = originalData[sortedPoints[start + leftRange / 2]][dim];
    auto thisNode = lsSmartPointer<treeNode>::New(level, dim, start, stop, median);
    treeNodes.push_back(thisNode);
    // thisNode->setRange(start, stop);
    // thisNode->level = level;
    // thisNode->dimSplit = dim;

    thisNode->left = buildByDepth(data, orders, originalData, start, start + leftRange, level + 1);
    thisNode->right = buildByDepth(data, orders, originalData, start + leftRange, stop, level + 1);
    //thisNode->isLeaf = (thisNode->left != nullptr || thisNode->right != nullptr ) ? false : true;

    return thisNode;
  }

  /*
  * ---- QUERY METHODS ---- *
  */

  /// [WIP] query the tree based on a single point
  ///
  /// checks if it should be a leaf treeNode
  /// TODO: TEST
  const lsSmartPointer<treeNode> nearest(const point_type &pt)
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning pt.")
          .print();
      return nullptr;
    }

    // TODO: Incorporate iterator
    lsSmartPointer<treeNode> thisNode = treeNodes.front();
    while (!thisNode->isLeaf)
      thisNode = (thisNode->belowMedian(pt)) ? thisNode->left : thisNode->right;
    return thisNode;
  }

  /*
  * ---- UTILITY METHODS ---- *
  */

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
    index_vector result(size);
    std::generate(result.begin(), result.end(), [n = 0]() mutable
                  { return n++; });
    // Taken and modified from: https://stackoverflow.com/questions/2312737/whats-a-good-way-of-temporarily-sorting-a-vector
    std::sort(result.begin(), result.end(),
              [toSortBegin, dim](auto left, auto right) {                // left, right : size_t -> indices of vector toSort
                return toSortBegin[left][dim] < toSortBegin[right][dim]; // sort in ascending order
              });
    return std::move(result);
  }

private:
  inline size_t getRoot(size_t treeNode)
  {
    return std::ceil((double)treeNode / 2) - (1);
  }

  inline size_t getLeftChild(size_t treeNode)
  {
    return treeNode * 2 + 1;
  }

  inline size_t getRightChild(size_t treeNode)
  {
    return treeNode * 2 + 2;
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
    for (size_t i = 0; i < treeNodes.size(); ++i)
    {
      std::string leaf = (treeNodes[i]->isLeaf) ? "#" : "";
      std::cout << "treeNode " << std::setw(3) << i << "(L" << treeNodes[i]->level << "): [" << treeNodes[i]->start << ", " << treeNodes[i]->stop << ")" << leaf;
      if (treeNodes[i]->isLeaf)
        std::cout << "#(" << treeNodes[i]->color << ")";
      std::cout << std::endl;
    }
  };

  /// prints tree one level per line at a time
  void printTreeByLevel()
  {
    size_t offset = 0;
    for (size_t level = 0; level <= depth; ++level)
    {
      size_t num_nodes = pow2(level);
      std::cout << "L" << level << ": ";
      for (size_t i = offset; i < offset + num_nodes; ++i)
      {
        //if (treeNodes[i]->level != level)
        //continue;
        if ((i + 1) % 8 == 0 && level > 3)
          std::cout << "\n";
        std::cout << std::setw(3) << i << "[" << treeNodes[i]->start << ", " << treeNodes[i]->stop << ") --- ";
      }
      std::cout << std::endl;
      offset += num_nodes;
    }
  };
};

#endif // LS_TREE
