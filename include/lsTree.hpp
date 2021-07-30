//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <cstdint>
#include <lsDimension.hpp>
#include <lsToDiskMesh.hpp>
#include <lsSmartPointer.hpp>

#define REORDER

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

constexpr int64_t divByPow2(int num, int exp)
{
  return num >> exp;
}

template <int D>
size_t inline convertToIndex(Dimension<D> dim)
{
  return dim.toArrayIndex();
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
  using Dim = Dimension<D>;
  using value_type = T;
  using size_type = std::size_t;

  using point_type = std::array<T, 3>; // always use point triple, even in 2D (D=2) case where z (3rd entry) is always zero
  using index_vector = std::vector<size_t>;
  using data_vector = std::vector<point_type> &; // TODO: Template this
  // ---------- treeNode ----------
  struct treeNode
  {
    bool isLeaf = true;
    size_t start = 0;
    size_t stop = 0;
    size_t level = 0;
    int color = 0;
    Dim dimSplit = Dim(D);
    std::string identifier = "";

    lsSmartPointer<treeNode> parent = nullptr;
    lsSmartPointer<treeNode> left = nullptr;
    lsSmartPointer<treeNode> right = nullptr;

    T median;

    treeNode()
    {
    }

    treeNode(size_t level_, size_t dim, size_t start_, size_t stop_, T median_, lsSmartPointer<treeNode> parent_) : level(level_), dimSplit(static_cast<Dim>(dim)), start(start_), stop(stop_), median(median_), parent(parent_)
    {
    }

    treeNode(size_t level_, size_t dim, size_t start_, size_t stop_, T median_) : level(level_), dimSplit(static_cast<Dim>(dim)), start(start_), stop(stop_), median(median_)
    {
    }

    treeNode(size_t level_, Dim dim, size_t start_, size_t stop_, T median_, lsSmartPointer<treeNode> parent_) : level(level_), dimSplit(dim), start(start_), stop(stop_), median(median_), parent(parent_)
    {
    }

    treeNode(size_t level_, Dim dim, size_t start_, size_t stop_, T median_) : level(level_), dimSplit(dim), start(start_), stop(stop_), median(median_)
    {
    }

    void setRange(size_t start_, size_t stop_)
    {
      start = start_;
      stop = stop_;
    }

    size_t size()
    {
      return stop - start;
    }

    size_t bucket_size()
    {
      return stop - start;
    }

    bool belowMedian(const point_type &pt)
    {
      return pt[dimSplit] <= median;
    }

    void assignColor(size_t colorAssign)
    {
      isLeaf = (left != nullptr || right != nullptr) ? false : true;
      color = colorAssign;
    }

    void identifyAs(std::string identifier_)
    {
      identifier = identifier_;
    }
  };

  using node_pointer = lsSmartPointer<treeNode>;
  using nodes_vector = std::vector<node_pointer>;
  using nodes_iterator = typename std::vector<node_pointer>::iterator;

  // ---------- DATA ----------
private:
  // ---------- Parameters + Attributes ----------
  std::string tree_type = "kd-Tree";
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  uint maxDepth = 4;
  uint maxNumBins = pow2(maxDepth);
  uint maxPointsPerBin = 4;
  uint depth = 0;
  uint numBins = 1; // numBins = 2^depth
  uint largestBinSize = 0;

public:
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  data_vector data;
  size_type N = 0; // number of points in data / size of data

  nodes_vector treeNodes;
  size_t startLeafs = 0; // pointer/iterator to start of 'leafs'
  std::vector<T> color;

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

    treeIterator(treeNode &node) : pos(&(this->data[node.start]))
    {
      auto bin_range = std::make_pair(&(this->data[node.start]), &(this->data[node.stop]));
      ranges = range_vector{bin_range};
    }
    treeIterator(nodes_vector &bins) : pos(&(this->data[bins[0].start]))
    {
      ranges = range_vector();
      for (const treeNode &node : bins)
      {
        auto bin_range = std::make_pair(&(this->data[node.start]), &(this->data[node.stop]));
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
      if (pos == ranges[segment].second)
      {
        ++segment;
        pos = ranges[segment].first;
        return *this;
      }
      ++pos;
      return *this;
    }

    treeIterator operator++(int)
    {
      auto temp = *this;
      if (temp.pos == temp.ranges[temp.segment][1])
      {
        ++temp.segment;
        temp.pos = temp.ranges[temp.segment][0];
        return temp;
      }
      ++temp.pos;
      return temp;
    }

    treeIterator &operator--()
    {
      if (pos == ranges[segment].first)
      {
        --segment;
        pos = ranges[segment].second;
        return &this;
      }
      --pos;
      return &this;
    }
    treeIterator operator--(int) //optional
    {
      auto temp = *this;
      if (temp.pos == temp.ranges[temp.segment].first)
      {
        --temp.segment;
        temp.pos = temp.ranges[temp.segment].second;
        return temp;
      }
      --temp.pos;
      return temp;
    }
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
#pragma endregion

#pragma region Interface
public:
  lsTree()
  {
  }

  lsTree(lsSmartPointer<lsMesh<T>> passedMesh)
      : mesh(passedMesh), data(passedMesh->getNodes()), N(passedMesh->getNodes().size()) {}

  void setMesh(lsSmartPointer<lsMesh<T>> passedMesh)
  {
    mesh = passedMesh;
    data = this->mesh->getNodes();
    N = data.size();
    treeNodes.clear();
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
    uint maxNumberOfNodes = pow2(maxDepth + 1) - 1;
    treeNodes.reserve(maxNumberOfNodes);
    // partition in x by median
    maxPointsPerBin = std::ceil((double)N / 16.);
    size_t medianPos = size_t(N / 2);

#ifndef NDEBUG // if in debug build
    {
      lsMessage::getInstance()
          .addDebug("lsTree: Building in DEBUG mode")
          .print();
      // lsVTKWriter<hrleCoordType>(surfaceMesh, lsFileFormatEnum::VTP,
      //                            "DEBUG_lsGeomAdvectMesh_contributewoMask.vtp")
      //     .apply();
      // auto mesh = lsSmartPointer<lsMesh<T>>::New();
      // lsToMesh<T, D>(maskLevelSet, mesh).apply();
      // lsVTKWriter<T>(mesh, lsFileFormatEnum::VTP,
      //                "DEBUG_lsGeomAdvectMesh_mask.vtp")
      //     .apply();
      // lsToMesh<T, D>(levelSet, mesh).apply();
      // lsVTKWriter<T>(mesh, lsFileFormatEnum::VTP,
      //                "DEBUG_lsGeomAdvectMesh_initial.vtp")
      //     .apply();
    }
#endif

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
      dimToSort++;
      auto y_idx = lsSmartPointer<index_vector>::New(argsortInDimension(data.begin(), dimToSort, N));
      orders.push_back(y_idx);
    }

    dimToSort++;
    auto z_idx = lsSmartPointer<index_vector>::New(N);
    std::generate(z_idx->begin(), z_idx->end(), [n = 0]() mutable
                  { return n++; });
    orders.push_back(z_idx);
    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2 -> y)
    sortedPoints = index_vector(N);
    std::generate(sortedPoints.begin(), sortedPoints.end(), [n = 0]() mutable
                  { return n++; });

    if (dimToSort != D - 1)
      std::cout << "OOPS" << std::endl;

    // Example
    // Global order
    // x: 0, 5, 3, 2, 1, 6, 4, 7,
    // y: 5, 1, 0, 4, 6, 3, 7, 2
    // z: 7, 3, 5, 2, 1, 4, 0, 6,
    // kd-Tree by level
    // L1(by x): sortedList = [0, 5,  3, 2] [1, 6,  4, 7]
    // L2(by y): sortedList = [5, 0] [3, 2] [1, 4] [6, 7]
    // L3(by z): sortedList = [5][0] [3][2] [1][4] [7][6]

    // std::cout << "D: " << D << std::endl;
    if (byLevel)
    {
      // Build the tree
      size_t binSize = N; // root node is parent to all data points
      // in-place construction of root
      treeNodes.emplace_back(lsSmartPointer<treeNode>::New(0, D, 0, N, data[sortedPoints[N / 2]][D]));

      for (size_t level = 1; (level < maxDepth + 1) && (binSize > maxPointsPerBin); ++level)
      {
        startLeafs = treeNodes.size();

        const size_t nodesAdded = buildLevel(sortedPoints, orders, level);

        // Compute largest bin of this Level
        for (auto nodeIt = beginLeafs(); nodeIt < treeNodes.end(); ++nodeIt)
        {
          binSize = (*nodeIt)->size() > binSize ? (*nodeIt)->size() : binSize;
        }
        largestBinSize = binSize;
      }
      numBins = pow2(depth);

      applyColorToNodes();

      sortByIndicesInplace(data, sortedPoints);

      // Assemble array of colors for each point
      color = std::vector<T>(N);
      for (auto it = beginLeafs(); it < treeNodes.end(); ++it)
      {
        auto treeNode = *it;
        for (size_t index = treeNode->start; index < treeNode->stop; ++index)
        {
          color[index] = treeNode->color;
        }
      }
      mesh->getPointData().insertNextScalarData(color, "color"); // adjustment for new ViennaLS
    }
    /// TODO: Implement and support byDepth build properly
    if (!byLevel)
    {
      auto root = buildByDepth(sortedPoints, orders, data, 0, N, 0);
    }
  }

private:
  // Dimension constexpr getNextDim(Dimension dim)
  // {
  //   Dimension next_dim(dim);
  //   ++next_dim;
  //   return next_dim;
  // }

  size_t constexpr getNextDimIdx(size_type level)
  {
    return D - 1 - level % (D);
  }

  // size_t constexpr getNextDimIdx(Dimension dim)
  // {
  //   Dimension next_dim = getNextDim(dim);
  //   return static_cast<size_t>(next_dim);
  // }

  /// Builds the tree level-by-level
  ///
  ///
  template <class Pointer>
  size_t buildLevel(index_vector &sortedPoints, std::vector<Pointer> &orders, size_t level)
  {
    depth = level; // new level --> set tree depth
    const size_t num_parents = treeNodes.size();
    const size_t startOfDirectParents = pow2(level - 1) - 1;

    // nodes_vector levelNodes;
    for (auto idx = startOfDirectParents; idx < num_parents; ++idx)
    {
      // root is to be split, so is not a leaf anymore
      lsSmartPointer<treeNode> root = treeNodes[idx];
      root->isLeaf = false;

      // Divide into two child nodes
      // Calculate start-, stop-indices and ranges
      const size_t start = root->start;
      const size_t stop = root->stop;
      const size_t range = stop - start;
      const size_t leftRange = range / 2;
      const size_t rightRange = range - leftRange;

      // Sort Range of the root
      //if (level > 1) // otherwise already sorted
      size_t sortDimIdx = convertToIndex(root->dimSplit);
      sortByDim(sortedPoints.begin() + start, sortedPoints.begin() + stop, sortDimIdx);

      Dim nextSplitDim = root->dimSplit--;
      // Make 2 new nodes
      // median = lookup Point from original data
      T median = data[sortedPoints[start + leftRange / 2]][sortDimIdx];
      root->left = lsSmartPointer<treeNode>::New(level, nextSplitDim, start, start + leftRange, median, root);

      median = data[sortedPoints[start + leftRange + rightRange / 2]][sortDimIdx];
      root->right = lsSmartPointer<treeNode>::New(level, nextSplitDim, start + leftRange, stop, median, root);

#ifndef NDEBUG // if in debug build
      {
        std::ostringstream leftId;
        leftId << root->dimSplit << " < " << root->median;
        root->left->identifyAs(leftId.str());
        std::ostringstream rightId;
        rightId << root->dimSplit << " > " << root->median;
        root->right->identifyAs(rightId.str());
      }
#endif

      // Insert pointers to new nodes into vector
      treeNodes.push_back(root->left);
      treeNodes.push_back(root->right);
    }
    // Returns however many nodes were added
    return treeNodes.size() - num_parents;
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

  //---- QUERY METHODS ----

private:
  lsSmartPointer<treeNode> getBin(const point_type &pt) const
  {
    lsSmartPointer<treeNode> thisNode = treeNodes.front();
    while (!thisNode->isLeaf)
      thisNode = (thisNode->belowMedian(pt)) ? thisNode->left : thisNode->right;
    return thisNode;
  }

public:
  node_pointer getRoot()
  {
    return treeNodes.front();
  }

  nodes_vector &getTreeNodes()
  {
    return treeNodes;
  }

  template <class iteratorType = iterator>
  iteratorType getNearestNeighbors(const point_type &point)
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iteratorType();
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
  iteratorType getNearestNeighbors(data_vector &points) const
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

  /*
  * ---- Configure TREE ---- *
  */

  void setMaxDepth(uint newDepth)
  {
    this->maxDepth = newDepth;
  }

  void setMaxNumBins(uint newLimit)
  {
    this->maxNumBins = newLimit;
  }

  void setMaxPointsPerBin(uint newLimit)
  {
    this->maxPointsPerBin = newLimit;
  }

  /// Reconfigure the trees parameters.
  ///
  /// Already built lsTree will not be rebuilt unless 'apply()' is called.
  void configureTreeSettings(uint maxDepth_, uint maxNumBins_, uint maxPointsPerBin_)
  {
    this->maxDepth = maxDepth_;
    this->maxNumBins = maxNumBins_;
    this->maxPointsPerBin = maxPointsPerBin_;
  }

  /*
  * ---- TREE INFO ---- *
  */
  size_type getNumberOfBins() const
  {
    return numBins;
  }

  size_type getNumberOfNodes() const
  {
    return treeNodes.size();
  }

  uint getDepth() const
  {
    return depth;
  }

  const std::string getTreeType()
  {
    return tree_type;
  }

  /// prints parameters
  void printInfo()
  {
    std::cout << getTreeType() << std::endl;
    std::cout << "lsTree Settings" << std::endl;
    std::cout << "depth: " << depth << std::endl;
    std::cout << "maxDepth: " << maxDepth << " --> depth: " << depth << std::endl;
    std::cout << "maxNumBins: " << maxNumBins << " --> numBins: " << numBins << std::endl;
    std::cout << "maxPointsPerBin: " << maxPointsPerBin << " --> largestBinSize: " << largestBinSize << std::endl;
  }

  void printBT(const std::string &prefix, lsSmartPointer<treeNode> node, bool isLeft)
  {
    if (node != nullptr)
    {
      std::cout << prefix;
      std::cout << (isLeft ? "├──" : "└──");

      if (node->identifier == "")
      {
        std::string dim = "z";
        if (node->dimSplit() == Dim::y)
          dim = "y";
        if (node->dimSplit() == Dim::x)
          dim = "x";
        std::cout << "(" << node->color << " / " << dim << " / " << node->median << ")" << std::endl;
      }
      else
      {
        std::cout << "(" << node->color << " / " << node->identifier << ")" << std::endl;
      }

      // enter the next tree level - left and right branch
      printBT(prefix + (isLeft ? "│   " : "    "), node->left, true);
      printBT(prefix + (isLeft ? "│   " : "    "), node->right, false);
    }
  }

  void printBT()
  {
    printBT("", treeNodes[0], false);
  }

#pragma region IteratorMethods
  iterator begin()
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iterator();
    }
    auto full_range = std::make_pair(&data.front(), &data.back());
    range_vector all{full_range};
    return iterator(all, &data.front());
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
    auto full_range = std::make_pair(&data.front(), &data.back());
    range_vector all{full_range};
    return iterator(all, &data.back());
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

private:
  nodes_iterator beginLeafs(void)
  {
    return std::next(treeNodes.begin(), startLeafs);
  }

  typename nodes_vector::iterator endLeafs(void)
  {
    return treeNodes.end();
  }
#pragma endregion All methods concerning treeIterator useage

#pragma Private Methods
private:
  /*
  * ---- UTILITY METHODS ---- *
  */

  void applyColorToNodes()
  {
    int b0 = numBins;
    auto &root = treeNodes.front();
    colorChild(root->left, root->color, b0, true);
    colorChild(root->right, root->color, b0, false);
  }

  void colorChild(lsSmartPointer<treeNode> node, size_t parentColor, int b0, bool isLeft)
  {
    int b = divByPow2(b0, node->level);
    node->color = isLeft ? parentColor - b : parentColor + b;
    if (node->left)
      colorChild(node->left, node->color, b0, true);
    if (node->right)
      colorChild(node->right, node->color, b0, false);
  };

  /// Sorts a vector (range) by the order given by second vector.
  /// The second vector contains the sorted list of indices of the first.
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
  /// order = [1,4,6,2]
  template <class It, class It2>
  void sortByIdxRange(It begin, It end, It2 beginSortBy, It2 endSortBy)
  {
    std::sort(begin, end,
              [beginSortBy, endSortBy](auto left, auto right)
              {
                // sort in ascending order
                for (auto iterator_to_sorted_data = beginSortBy; iterator_to_sorted_data < endSortBy; ++iterator_to_sorted_data)
                {
                  if (*iterator_to_sorted_data == left)
                    return true; // switch order of left and right
                  if (*iterator_to_sorted_data == right)
                    return false; // left and right are already ordered correctly
                }
                return false; // neither could be found...should this give an error?
              });
  }
  template <class It>
  void sortByDim(It begin, It end, uint dim)
  {
    std::sort(begin, end,
              [this, &dim](const auto &left, const auto &right)
              {
                return data[left][dim] < data[right][dim]; // false...switch order of left and right
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
