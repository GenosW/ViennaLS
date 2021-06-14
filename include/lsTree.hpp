//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <cstdint>
#include <lsToDiskMesh.hpp>
#include <lsSmartPointer.hpp>

#define REODER

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
/// TODO: order
///
/// Types:
///
/// using value_type = T;
/// using size_type = std::size_t
/// using point_type = std::array<T, 3>;
/// using index_vector = std::vector<size_t>;
/// using data_vector = std::vector<point_type>;
/// using iterator = treeIterator<point_type>;
/// using const_iterator = constTreeIterator<point_type>;
/// using range_vector = typename iterator::range_vector;
/// using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
template <class T, int D>
class lsTree
{
  // ---------- Typedefs etc. ----------
public:
  using value_type = T;
  using size_type = std::size_t;

  using point_type = std::array<T, 3>;
  using index_vector = std::vector<size_t>;
  using data_vector = std::vector<point_type>; // TODO: Template this

private:
  // ---------- Parameters + Attributes ----------
  std::string tree_type = "kd-Tree";
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  uint maxDepth = 4;
  uint maxNumBins = pow2(maxDepth);
  uint maxPointsPerBin = 4;

  // ---------- treeNode ----------
#pragma region treeNode
  struct treeNode
  {
    bool isLeaf = false;
    /// TODO: Implement std::pair<pointer, pointer> for start stop?
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
    /// same as size()
    ///
    /// TODO: stdlib bucket type containers interface?
    size_t bucket_size()
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

  // ---------- DATA ----------
private:
  using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
  // The mesh we're building a tree for
  lsSmartPointer<lsMesh<T>> mesh = nullptr;

  uint depth = 0;
  uint numBins = 1; // numBins = 2^depth

  nodes_vector treeNodes;
  data_vector data;
  std::vector<T> color;
#pragma endregion

#pragma region treeIterator
  /// Pointer based implementation
  ///
  /// Nicely compiled list of needed methods:
  /// https://stackoverflow.com/questions/7758580/writing-your-own-stl-container/7759622#7759622
  ///
  /// Some points about vector iterators:
  /// https://codereview.stackexchange.com/questions/202157/basic-iterator-supporting-vector-implementation
  /// Iterators: https://www.cplusplus.com/reference/iterator/
  template <class VT = point_type>
  struct treeIterator
  {
  public:
    // TODO: Evaluate iterator category
    using iterator_category = std::bidirectional_iterator_tag; // seems likeliest for now
    // std::random_access_iterator_tag;
    // std::forward_iterator_tag;
    using value_type = VT;
    using difference_type = size_type;
    using pointer = value_type *;
    using reference = value_type &;

  public:
    using range_vector = const std::vector<std::pair<pointer, pointer>>;

    // Position based implementaiton
    range_vector ranges;
    size_t segment = 0;
    // Idea: Utilize [start, stop) of each leaf node
    // ranges = [(startBin1, stopBin1), (startBin2, stopBin2),
    //           (startBin3, stopBin3), (startBin4, stopBin4), ... ]
    pointer pos; // ->data[init]
    // unordered map/set iterator anschauen

  public:
    treeIterator() = default;

    treeIterator(range_vector input) : pos(input[0].first), ranges(input){};
    treeIterator(range_vector input, size_type start) : pos(input[start].first), ranges(input){};
    treeIterator(range_vector input, pointer start_pos) : pos(start_pos), ranges(input){};
    treeIterator(treeIterator &other) : pos(other.pos), ranges(other.ranges){};

    treeIterator(treeNode &node) : pos(*(this->data[node.start]))
    {
      auto bin_range = std::make_pair(*(this->data[node.start]), *(this->data[node.stop]));
      ranges = range_vector{bin_range};
    }
    treeIterator(nodes_vector &bins) : pos(*(this->data[bins[0].start]))
    {
      ranges = range_vector();
      for (const treeNode &node : bins)
      {
        auto bin_range = std::make_pair(*(this->data[node.start]), *(this->data[node.stop]));
        ranges.push_back(bin_range);
      }
    }

    // treeIterator &operator=(treeIterator &other)
    // {
    //   return treeIterator(other);
    // }

    // ------- Access -------
    reference operator*() const
    {
      return *pos;
    }

    pointer operator->() const
    {
      return &(*pos);
    }

    //reference operator[](size_type) const //optional

    // ---------- Comparison Operations ----------
    bool operator==(const treeIterator &other) const
    {
      return this->pos == other.pos;
    }

    bool operator!=(const treeIterator &other) const
    {
      return this->pos != other.pos;
    }
    bool operator<(const treeIterator &other) const //optional
    {
      return this->pos < other.pos;
    }
    bool operator>(const treeIterator &other) const //optional
    {
      return this->pos > other.pos;
    }
    bool operator<=(const treeIterator &other) const //optional
    {
      return this->pos <= other.pos;
    }
    bool operator>=(const treeIterator &other) const //optional
    {
      return this->pos >= other.pos;
    }

    // ---------- Arithmetic Operations ----------
    treeIterator &operator++()
    {
      if (pos == ranges[segment][1])
      {
        ++segment;
        pos = ranges[segment][0];
        return &this;
      }
      ++pos;
      return &this;
    }

    treeIterator operator++(int)
    {
      if (pos == ranges[segment][1])
      {
        ++segment;
        pos = ranges[segment][0];
        return &this;
      }
      pos++;
      return &this;
    }

    treeIterator &operator--()
    {
      if (pos == ranges[segment][0])
      {
        --segment;
        pos = ranges[segment][1];
        return &this;
      }
      --pos;
      return &this;
    }
    treeIterator operator--(int) //optional
    {
      if (pos == ranges[segment][0])
      {
        --segment;
        pos = ranges[segment][1];
        return &this;
      }
      pos--;
      return &this;
    }

    /// TODO: depends on iterator type...
    /// These types of operations might not be possible due to disconnected segments
    // treeIterator &operator+=(const size_type n) //optional
    // {
    //   this->pos += n;
    //   return &this;
    // }
    // friend treeIterator operator+(const size_type n, const treeIterator &); //optional

    // treeIterator &operator-=(const size_type n) //optional
    // {
    //   this->pos -= n;
    //   return &this;
    // }
    // TODO: Evaluate if needed/usefull/sensible
    // treeIterator operator+(const size_type n) const //optional
    // {
    //   treeIterator ret(this);
    //   ret.pos += n;
    //   return ret;
    // }
    // treeIterator operator-(const size_type n) const //optional
    // {
    //   treeIterator ret(this);
    //   ret.pos -= n;
    //   return ret;
    // }

    /// TODO: Difference could still be useful, but maybe only within segments?
    // difference_type operator-(const size_type n) const //optional
    // {
    //   this->pos + n;
    //   return this;
    // }
  };

  /// const_iterator version of the (tree)Iterator
  template <class VT = point_type>
  struct constTreeIterator : treeIterator<const VT>
  {
  public:
    // MyConstIterator() = default; // default constructor enabled
    //using treeIterator<const VT>::treeIterator;

    /// Construct from non-const iterator
    constTreeIterator(const treeIterator<VT> &nonconst_other) : treeIterator<VT>(nonconst_other){};
  };

  using iterator = treeIterator<point_type>;
  using const_iterator = constTreeIterator<point_type>;
  using range_vector = typename iterator::range_vector;
// treeIterator() = default;
// treeIterator(size_type start) : pos(*data[start]){};
// treeIterator(treeNode *init) : pos(*(init->start)){};
// treeIterator(range_vector input) : pos(input[0].first), ranges(input){};
// treeIterator(range_vector input, size_type start) : pos(input[start].first), ranges(input){};
// treeIterator(range_vector input, pointer start_pos) : pos(start_pos), ranges(input){};
#pragma endregion

#pragma region Interface
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

    auto data = mesh->getNodes();
    index_vector x;
    index_vector y;
    index_vector z;
    index_vector sortedPoints;

    std::vector<lsSmartPointer<index_vector>> orders;

    // 1.Generate absolute orders in each dimension
    // Sort in z/y axis (D=3/D=2)
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    // Sort in x dimension (always)
    size_t dimToSort = 0;
    auto x_idx = lsSmartPointer<index_vector>::New(argsortInDimension(data.begin(), dimToSort, N));
    orders.push_back(x_idx);

    // Sort in y dimension (if D = 3)
    if constexpr (D > 2)
    {
      dimToSort++; // dimToSort = 1 <==> y
      auto y_idx = lsSmartPointer<index_vector>::New(largsortInDimension(data.begin(), dimToSort, N));
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
      // Build the tree
      size_t binSize = N;
      auto root = lsSmartPointer<treeNode>::New(0, dimToSort, 0, N, data[sortedPoints[N / 2]][dimToSort]);
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
        startLastLevel = buildLevel(sortedPoints, orders, data, level, dim, N);
        binSize = treeNodes.back()->size();
      }
      numBins = pow2(depth);
      // Color the data nodes
      size_t idx = 0;
      // for (auto &treeNode : treeNodes)
      // {
      //   treeNode->color = idx;
      //   ++idx;
      // }

      for (size_t i = startLastLevel; i < treeNodes.size(); ++i)
      {
        treeNodes[i]->color = idx;
        ++idx;
      }
#define REORDER
#ifdef REORDER
      //data = [(x,y,z), (x,y,z), (x,y,z),]

      //sortedPoints = [0, 1, 2, 90,91,92 ||, ...]
      // TODO: impact on large vs small data arrays
      sortByIndicesInplace(data, sortedPoints);
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

#endif
#ifndef REORDER
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
#endif
      mesh->insertNextScalarData(color, "color");
    }
    /// TODO: Implement and support byDepth build properly
    if (!byLevel)
    {
      auto root = buildByDepth(sortedPoints, orders, data, 0, N, 0);
    }
  }

  //---- QUERY METHODS ----

private:
  /// Helper function
  lsSmartPointer<treeNode> getBin(const point_type &pt) const
  {

    lsSmartPointer<treeNode> thisNode = treeNodes.front();
    while (!thisNode->isLeaf)
      thisNode = (thisNode->belowMedian(pt)) ? thisNode->left : thisNode->right;
    return thisNode;
  }

public:
  iterator getNearestNeighbors(const point_type &point)
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iterator();
    }
    return iteratorType(getBin(point));
  }

  /// Get the nearest Neighbors in form of one or more bins.
  ///
  /// Returns an iterator to the nearest points found.
  /// Since the lsTree partioned the space into multiple bins, the resulting
  /// iterator of this query might refer to multiple disconnected segments
  /// in memory - even if the original container was contingent
  /// (e.g. std::vector, std::array).
  template <class iteratorType = iterator>
  iteratorType getNearestNeighbors(const data_vector &points)
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iteratorType();
    }

    nodes_vector found_nodes;
    for (auto &pt : points)
      found_nodes.push_back(getBin(pt));

    return iteratorType(found_nodes);
  }
#pragma region IteratorMethods
  // private:
  // template <class returnType>
  // returnType treeBuiltCheck(returnType& toReturnIfError){
  //   if (mesh == nullptr)
  //     {
  //       lsMessage::getInstance()
  //           .addWarning("No mesh was passed to lsTree.")
  //           .print();
  //       return toReturnIfError;
  //     }
  // }

  iterator begin()
  {

    range_vector all{data.front(), data.back()};
    return iterator(all);
  };
  iterator end()
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iterator();
    }
    range_vector all{data.front(), data.back()};
    return const_iterator(all, data.back());
  }
  const_iterator begin() const
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return const_iterator();
    }
    range_vector all{data.front(), data.back()};
    return const_iterator(all);
  };
  const_iterator end() const
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return const_iterator();
    }
    range_vector all{data.front(), data.back()};
    return const_iterator(all, data.back());
  };

  const_iterator cbegin() const
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return const_iterator();
    }
    range_vector all{data.front(), data.back()};
    return const_iterator(all);
  };
  const_iterator cend() const
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return const_iterator();
    }
    range_vector all{data.front(), data.back()};
    return const_iterator(all, data.back());
  };
#pragma endregion All methods concerning treeIterator useage
  /*
  * ---- TREE INFO ---- *
  */

  size_type getNumberOfBins() const
  {
    return numBins;
  }
  /// stdlib equivalent
  ///
  /// TODO: stdlib bucket type containers interface?
  size_type bucket_count() const
  {
    return numBins;
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
  }

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
  }

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
  }
#pragma endregion These methods define the interactive(method) interface to the outside.
#pragma Private Methods
private:
  /*
  * ---- UTILITY METHODS ---- *
  */

  /// Builds the tree level-by-level
  ///
  ///
  /// TODO: Use Iterators (begin/end) instead of start/stop
  template <class Vector, class VectorOfPointers, class VectorData>
  size_t buildLevel(Vector &data, VectorOfPointers &orders, VectorData &originalData, size_t level, size_t dim, size_t N)
  {
    size_t num_nodes = pow2(level);

    depth = level; // new level --> set tree depth
    std::cout << "Level: " << level << " [dim=" << dim << "]" << std::endl;
    std::cout << "num_nodes: " << num_nodes << std::endl;
    std::cout << "nodeSize: ca. " << (N / num_nodes) << std::endl;
    std::cout << "treeNodes.size(): " << treeNodes.size() << std::endl;
    std::cout << "parent nodes: " << treeNodes.size() - (pow2(level - 1) - 1) << std::endl;

    auto order_begin = orders[dim]->begin();
    auto order_end = orders[dim]->end();
    auto data_begin = data.begin();
    auto data_end = data.end();

    size_t num_parents = treeNodes.size();

    nodes_vector levelNodes;
    // for (auto rootsIt = treeNodes.begin() + pow2(level - 1) - 1; rootsIt < treeNodes.end(); ++rootsIt)
    for (auto idx = pow2(level - 1) - 1; idx < num_parents; ++idx)
    {
      std::cout << "Root " << idx << " - ";
      // root is to be split, so is not a leaf anymore
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
      // median = lookup Point from original data
      T median = originalData[data[start + leftRange / 2]][dim];
      lsSmartPointer<treeNode> left = lsSmartPointer<treeNode>::New(level, dim, start, start + leftRange, median);
      root->left = left;

      median = originalData[data[start + leftRange + rightRange / 2]][dim];
      lsSmartPointer<treeNode> right = lsSmartPointer<treeNode>::New(level, dim, start + leftRange, stop, median);
      root->right = right;

      // Insert pointers to new nodes into vector
      levelNodes.push_back(left);
      levelNodes.push_back(right);
    }

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

    T median = originalData[data[start + leftRange / 2]][dim];
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

  /// Builds a new vector out of the original data and desired sorting order (index array)
  template <class Vector, class VectorIndices>
  Vector sortByIndices(Vector &data, VectorIndices &indices)
  {
    if (data.size() != indices.size())
      std::cout << "sizes do not match!" << std::endl;

    size_t size = data.size();
    Vector result(size);
    for (int i = 0; i < size; ++i)
    {
      result[i] = data[indices[i]];
    }
    return std::move(result);
  }
  /// Builds a new vector out of the original data and desired sorting order (index array)
  // template <class Vector, class VectorIndices>
  // Vector sortByIndices2(Vector &data, VectorIndices &indices)
  // {
  //   using vt = Vector::value_type;
  //   if (data.size() != indices.size())
  //     std::cout << "sizes do not match!" << std::endl;

  //   size_t size = data.size();
  //   bool tempHot = false;
  //   vt temp;
  //   for (int i = 0; i < size; ++i)
  //   {
  //     temp = ()
  //         data[i] = data[indices[i]];
  //   }
  //   return std::move(result);
  // }

  /// Sorts a vector data in place by ordering given via vector indices (same as sortByIndices but in place)
  template <class Vector, class VectorIndices>
  void sortByIndicesInplace(Vector &data, VectorIndices &indices)
  {
    if (data.size() != indices.size())
      std::cout << "sizes do not match!" << std::endl;

    size_t size = data.size();
    Vector result(size);
    for (int i = 0; i < size; ++i)
    {
      result[i] = data[indices[i]];
    }
    data = std::move(result);
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
};

#endif // LS_TREE
