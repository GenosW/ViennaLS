//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <numeric>
#include <lsDimension.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>
#include <vector>

// #define REORDER

namespace lsInternal
{
  constexpr int64_t ipow(int64_t base, int exp, int64_t result = 1)
  {
    return exp < 1
               ? result
               : ipow(base * base, exp / 2, (exp % 2) ? result * base : result);
  }

  // all depth values are known at compile time
  constexpr int64_t pow2(int64_t exp) { return 1 << exp; }

  // fast (constexpr) function to calculate: num / 2^exp
  constexpr int64_t divByPow2(int64_t num, int64_t exp) { return num >> exp; }

  template <int D>
  std::size_t inline convertToIndex(lsInternal::Dimension<D> dim)
  {
    return dim.toArrayIndex();
  }

  template <typename Point>
  void print_point(Point &p, std::ostream &os = std::cout)
  {
    os << "[" << p[0] << ", " << p[1] << ", " << p[2] << "]";
    // so << "(" << *p[0][0] << "," << *p[0][1] << "," << *p[0][2] << ")" << std::endl;
  };

#ifndef __cpp_lib_interpolate
  /// Implementation of C++20 std::midpoint (defined in <numeric>)
  /// in case it use not implemented (e.g. when compiling with C++17)
  template <class T>
  constexpr std::enable_if_t<std::is_floating_point_v<T>, T>
  midpoint(T a, T b) noexcept
  {
    if (std::isnan(a) || std::isnan(b))
      return std::numeric_limits<T>::quiet_NaN();
    else
      return std::isnormal(a) && std::isnormal(b)
                 ? a / 2 + b / 2
                 : (a + b) / 2;
  };
#else
  using std::midpoint;
#endif
} // end namespace lsInternal

/// Tree structure that seperates a domain
///
///
/// Types:
///
/// using Dim = typename lsInternal::Dimension<D>;
/// using value_type = T;
/// using size_type = std::size_t;
/// using point_type = std::array<T, 3>;
/// using Bin = std::vector<point_type>;
/// using Data_Container = std::vector<Bin>;
/// using iterator = typename Bin::iterator;
/// using const_iterator = typename Bin::const_iterator;
/// using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
/// using nodes_iterator = typename std::vector<lsSmartPointer<treeNode>>::iterator;
/// using nodes_const_iterator = typename std::vector<lsSmartPointer<treeNode>>::const_iterator;

template <class T, int D>
class lsTree
{
  // ---------- Typedefs etc. ----------
public:
  using Dim = typename lsInternal::Dimension<D>;
  using value_type = T;
  using size_type = std::size_t;

  using point_type = std::array<T, 3>;
  // always use point triple, even in 2D (D=2) case where
  // z (3rd entry) is always zero
  using Bin = std::vector<point_type>;
  using Data_Container = std::vector<Bin>; // TODO: Template this
  using iterator = typename Bin::iterator;
  using const_iterator = typename Bin::const_iterator;

  // ---------- treeNode ----------
  struct treeNode
  {
    // lsTree &owner;
    bool isLeaf = true;
    size_type start = 0;
    size_type stop = 0;
    size_type level = 0;
    size_type leafNum = 0;
    int color = 0;
    Dim dimensionToSplit = Dim(D);
    std::string identifier = "";

    lsSmartPointer<treeNode> parent = nullptr;
    lsSmartPointer<treeNode> left = nullptr;
    lsSmartPointer<treeNode> right = nullptr;

    T median;

    treeNode() = default;

    treeNode(size_type level_, size_type dim, size_type start_, size_type stop_, T median_,
             lsSmartPointer<treeNode> parent_)
        : level(level_), dimensionToSplit(static_cast<Dim>(dim)), start(start_),
          stop(stop_), median(median_), parent(parent_) {}

    treeNode(size_type level_, size_type dim, size_type start_, size_type stop_, T median_)
        : level(level_), dimensionToSplit(static_cast<Dim>(dim)), start(start_),
          stop(stop_), median(median_) {}

    treeNode(size_type level_, Dim dim, size_type start_, size_type stop_, T median_,
             lsSmartPointer<treeNode> parent_)
        : level(level_), dimensionToSplit(dim), start(start_), stop(stop_),
          median(median_), parent(parent_) {}

    treeNode(size_type level_, Dim dim, size_type start_, size_type stop_, T median_)
        : level(level_), dimensionToSplit(dim), start(start_), stop(stop_),
          median(median_) {}

    void setRange(size_type start_, size_type stop_)
    {
      start = start_;
      stop = stop_;
    }

    size_type size() { return stop - start; }

    size_type bucket_size() { return stop - start; }

    bool belowMedian(const point_type &pt)
    {
      return pt[dimensionToSplit] <= median;
    }

    void assignColor(size_type colorAssign)
    {
      isLeaf = (left != nullptr || right != nullptr) ? false : true;
      color = colorAssign;
    }

    void identifyAs(std::string identifier_) { identifier = identifier_; }

    friend std::ostream &operator<<(std::ostream &os, const treeNode &node)
    {
      os << node.color << ": [" << node.start << ", " << node.top << "]";
      return os;
    }

    iterator begin(lsTree &tree)
    {
      return tree.data[leafNum].begin();
    }

    iterator end(lsTree &tree)
    {
      return tree.data[leafNum].end();
    }
  };

  using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
  using nodes_iterator = typename std::vector<lsSmartPointer<treeNode>>::iterator;
  using nodes_const_iterator = typename std::vector<lsSmartPointer<treeNode>>::const_iterator;

private:
  // ---------- Parameters ----------
  std::string tree_type = "kd-Tree";
  bool frozen = false;
  bool dataGiven = false;
  bool useBestSplit = false;
  size_type maxDepth = 4; // TODO: maxDepth prop log10(N)
  size_type maxNumBins = lsInternal::pow2(maxDepth);
  size_type maxPointsPerBin = 4;

public:
  Data_Container data;
  // bins == std::vector<point>
  size_type depth = 0;
  size_type N = 0;       // number of points in data / size of data
  size_type numBins = 1; // numBins = 2^depth
  size_type largestBinSize = 0;

  lsSmartPointer<treeNode> root = nullptr;
  nodes_vector treeNodes;
  size_type startLeafs = 0; // TODO: Make pointer/iterator to start of 'leafs' instead? --> consider consequences!
  std::vector<T> color;

  // TODO: Figure out which typing to use here...
  // using iterator = typename Data_Container::iterator;
  // using const_iterator = typename Data_Container::const_iterator;

#pragma region IteratorMethods
  iterator begin()
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iterator(*this);
    }
    return data.begin(); // FIXME: just to show what should come out
  };
  iterator end()
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iterator(*this);
    }
    return data.end(); // FIXME: just to show what should come out
  }

  // The following iterator (interfaces) aren't yet implemented.
  // These just represent a general look. 

  // const_iterator begin() const
  // {
  //   if (treeNodes.empty())
  //   {
  //     lsMessage::getInstance()
  //         .addWarning("lsTree was not built. Returning empty iterator.")
  //         .print();
  //     return const_iterator();
  //   }
  //   range_vector all{data.front(), data.back()};
  //   return const_iterator(all);
  // };
  // const_iterator end() const
  // {
  //   if (treeNodes.empty())
  //   {
  //     lsMessage::getInstance()
  //         .addWarning("lsTree was not built. Returning empty iterator.")
  //         .print();
  //     return const_iterator();
  //   }
  //   range_vector all{data.front(), data.back()};
  //   return const_iterator(all, data.back());
  // };

  // const_iterator cbegin() const
  // {
  //   if (treeNodes.empty())
  //   {
  //     lsMessage::getInstance()
  //         .addWarning("lsTree was not built. Returning empty iterator.")
  //         .print();
  //     return const_iterator();
  //   }
  //   range_vector all{data.front(), data.back()};
  //   return const_iterator(all);
  // };
  // const_iterator cend() const
  // {
  //   if (treeNodes.empty())
  //   {
  //     lsMessage::getInstance()
  //         .addWarning("lsTree was not built. Returning empty iterator.")
  //         .print();
  //     return const_iterator();
  //   }
  //   range_vector all{data.front(), data.back()};
  //   return const_iterator(all, data.back());
  // };

private:
  // Some helper methods to access only the leaf nodes

  nodes_iterator beginLeafs(void)
  {
    return std::next(treeNodes.begin(), startLeafs);
  }

  nodes_const_iterator beginLeafs(void) const
  {
    return std::next(treeNodes.begin(), startLeafs);
  }

  nodes_iterator endLeafs(void) { return treeNodes.end(); }

  nodes_const_iterator endLeafs(void) const { return treeNodes.end(); }

#pragma endregion All methods concerning treeIterator useage

public:
  lsTree() {}

  lsTree(lsSmartPointer<lsMesh<T>> passedMesh) : N(passedMesh->getNodes().size())
  {

    data = std::vector<Bin>{passedMesh->getNodes()};
    dataGiven = true;
#ifndef NDEBUG // if in debug build
    {
      lsMessage::getInstance()
          .addDebug("lsTree: used Constructor lsTree(lsSmartPointer<lsMesh<T>> passedMesh)")
          .print();

      lsMessage::getInstance()
          .addDebug("lsTree: data[0][0][0] = " + std::to_string(data[0][0][0]))
          .print();
    }
#endif
  }

  lsTree(lsMesh<T> const &passedMesh)
      : N(passedMesh.getNodes().size())
  {
    data = std::vector<Bin>{passedMesh.getNodes()};
    dataGiven = true;
#ifndef NDEBUG // if in debug build
    {
      lsMessage::getInstance()
          .addDebug("lsTree: used Constructor lsTree(lsMesh<T> const &passedMesh)")
          .print();

      lsMessage::getInstance()
          .addDebug("lsTree: data[0][0][0] = " + std::to_string(data[0][0][0]))
          .print();
    }
#endif
  }

  /// entry point
  lsTree &apply(bool byLevel = true)
  {
    if (!dataGiven)
    {
      lsMessage::getInstance()
          .addWarning("No data was passed to lsTree! Did not do anything.")
          .print();
      return *this;
    }
    uint maxNumberOfNodes = lsInternal::pow2(maxDepth + 1) - 1;
    treeNodes.reserve(maxNumberOfNodes);
    // partition in x by median
    maxPointsPerBin = std::ceil((double)N / 16.);
    size_type medianPos = size_type(N / 2);

#ifndef NDEBUG // if in debug build
    {
      lsMessage::getInstance()
          .addDebug("lsTree: Building in DEBUG mode")
          .print();
    }
#endif

    std::vector<size_type> sortedPoints;

    // Nodes are already sorted correctly in highest dimension (D=3 -> z / D=2
    // -> y)
    sortedPoints = std::vector<size_type>(N);
    std::generate(sortedPoints.begin(), sortedPoints.end(),
                  [n = 0]() mutable
                  { return n++; });

    // Example
    // Global order
    // x: 0, 5, 3, 2, 1, 6, 4, 7,
    // y: 5, 1, 0, 4, 6, 3, 7, 2
    // z: 7, 3, 5, 2, 1, 4, 0, 6,
    // kd-Tree by level
    // L1(by x): sortedList = [0, 5,  3, 2] [1, 6,  4, 7]
    // L2(by y): sortedList = [5, 0] [3, 2] [1, 4] [6, 7]
    // L3(by z): sortedList = [5][0] [3][2] [1][4] [7][6]

    if (byLevel)
    {
      // Build the tree
      auto binSize = N;
      // in-place construction of root
      // "highest" coord = D - 1
      // e.g. D = 2 --> coord = 1 ^= y
      // e.g. D = 3 --> coord = 2 ^= z
      root = lsSmartPointer<treeNode>::New(0, Dim(D - 1), 0, N, data[0][sortedPoints[N / 2]][D - 1]);

#ifndef NDEBUG // if in debug build
      {
        root->identifyAs("root");
      }
#endif
      // treeNodes.emplace_back(root);
      treeNodes.push_back(root);

      for (size_type level = 1; (level < maxDepth + 1) && (binSize > maxPointsPerBin); ++level)
      {
        startLeafs = treeNodes.size(); // set to end when new level is built

        const auto nodesAdded = buildLevel(sortedPoints, level);

        // Compute largest bin of this Level
        auto maxBin = 0;
        for (auto nodeIt = beginLeafs(); nodeIt < endLeafs(); ++nodeIt)
          maxBin = (*nodeIt)->size() > maxBin ? (*nodeIt)->size() : maxBin;
        largestBinSize = maxBin;
      }

      numBins = lsInternal::pow2(depth);
#ifndef NDEBUG // if in debug build
      {
        applyColorToNodes(); // Assemble array of colors for each point
        // TODO: Refactor assembly of color
        color = std::vector<T>(N);
        for (auto it = beginLeafs(); it < treeNodes.end(); ++it)
        {
          auto treeNode = *it;
          for (size_type idx = treeNode->start; idx < treeNode->stop; ++idx)
          {
            unsigned originalIndex = sortedPoints[idx];
            color[originalIndex] = treeNode->color;
          }
        }
      }
#endif
      sortByIndicesInplace(data, sortedPoints);
    }
    /// TODO: Implement and support byDepth build properly
    if (!byLevel)
    {
      // auto root = buildByDepth(sortedPoints, data[0], 0, N, 0);
    }

    // We freeze the tree once it is done building.
    // This tree does not support insertion of additional data,
    // since we assume that it anyway needs to be rebuilt for every change
    // in the underlying data. I.e. we assume "massive" changes.
    frozen = true;
    return *this;
  }

private:
  /// Builds the tree level-by-level
  ///
  ///
  size_type buildLevel(std::vector<size_type> &sortedPoints, size_type level)
  {
    depth = level; // new level --> set tree depth
    const auto num_parents = treeNodes.size();
    const auto startOfDirectParents = lsInternal::pow2(level - 1) - 1;

    // nodes_vector levelNodes;
    for (auto idx = startOfDirectParents; idx < num_parents; ++idx)
    {
      // root is to be split, so is not a leaf anymore
      lsSmartPointer<treeNode> root = treeNodes[idx];
      root->isLeaf = false;

      // Divide into two child nodes
      // Calculate start-, stop-indices and ranges
      const auto start = root->start;
      const auto stop = root->stop;
      const auto range = stop - start;
      const auto leftRange = range / 2;
      const auto rightRange = range - leftRange;

      // Figure where the sorting should take place!

      // Sort Range of the root
      // if (level > 1) // otherwise already sorted
      Dim nextSplitDim;
      if (this->useBestSplit)
      {
        unsigned nextIdx = computeNextBestSplitDimension(sortedPoints.begin() + start, sortedPoints.begin() + stop);
        Dim nextSplitDim(nextIdx);
      }
      else
      {
        nextSplitDim = root->dimensionToSplit--;
      }
      unsigned sortDimIdx = lsInternal::convertToIndex(nextSplitDim);

      // Make 2 new nodes
      unsigned leftMedianIndex = start + leftRange / 2;
      unsigned rightMedianIndex = start + leftRange + rightRange / 2;

      T leftAverage = 0, rightAverage = 0;

      // [start:stop] --> [start : leftStop] | [rightstart : stop]
      //              --> [start : leftStop] | [leftStop   : stop]
      const auto leftStart = start, leftStop = start + leftRange;
      const auto rightStart = leftStop, rightStop = stop;
      sortByDim(sortedPoints.begin() + leftStart, sortedPoints.begin() + leftStop, sortDimIdx);
      sortByDim(sortedPoints.begin() + rightStart, sortedPoints.begin() + rightStop, sortDimIdx);

      // median = lookup Point from original data
      auto get_median = [this, &sortedPoints, sortDimIdx](size_type medianIndex, size_type range)
      {
        // auto result = (range % 2 == 0) ? (
        //                                      (data[0][sortedPoints[medianIndex]][sortDimIdx] +
        //                                       data[0][sortedPoints[medianIndex - 1]][sortDimIdx]) /
        //                                      2)
        //                                : (data[0][sortedPoints[medianIndex]][sortDimIdx]);
        auto result = (range % 2 == 0) ? lsInternal::midpoint(data[0][sortedPoints[medianIndex]][sortDimIdx], data[0][sortedPoints[medianIndex - 1]][sortDimIdx])
                                       : (data[0][sortedPoints[medianIndex]][sortDimIdx]);
        return result;
      };

      leftAverage = get_median(leftMedianIndex, leftRange);
      rightAverage = get_median(rightMedianIndex, rightRange);

      root->left = lsSmartPointer<treeNode>::New(level, nextSplitDim, leftStart, leftStop, leftAverage, root);
      root->right = lsSmartPointer<treeNode>::New(level, nextSplitDim, rightStart, rightStop, rightAverage, root);

#ifndef NDEBUG // if in debug build
      {
        std::ostringstream leftId;
        leftId << root->dimensionToSplit << " < " << root->median;
        root->left->identifyAs(leftId.str());
        std::ostringstream rightId;
        rightId << root->dimensionToSplit << " > " << root->median;
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

  //---- QUERY METHODS ----

  // private: // TODO: make it private again
public:
  /// method to get the node, that represents the neighborhood bin of a point.
  lsSmartPointer<treeNode> getBin(const point_type &pt, bool trace = false) const
  {
    lsSmartPointer<treeNode> thisNode = treeNodes.front();
    while (!thisNode->isLeaf)
    {
      if (trace)
      {
        std::stringstream message;
        message << "node(color=" << thisNode->color << ", median: " << thisNode->dimensionToSplit << "=" << thisNode->median << ", id: " << thisNode->identifier << ")";
        lsMessage::getInstance().addDebug(message.str()).print();
      }
      thisNode = (thisNode->belowMedian(pt)) ? thisNode->left : thisNode->right;
    }
    return thisNode;
  }

  template <typename Tv, int Dv>
  lsSmartPointer<treeNode> getBin(const hrleVectorType<Tv, Dv> &vector, bool trace = false) const
  {
    if constexpr (Dv < 3)
    {
      point_type pt{vector[0], vector[1], 0};
      return this->getBin(pt, trace);
    }

    point_type pt{vector[0], vector[1], vector[2]};
    return this->getBin(pt, trace);
  }

public:
  bool isFrozen()
  {
    return frozen;
  }

  template <typename T_Mesh>
  void addColor(lsSmartPointer<lsMesh<T_Mesh>> mesh)
  {
    mesh->getPointData().insertNextScalarData(color, "color");
    // adjustment for new ViennaLS
  }

  template <typename T_Mesh>
  void addColor(lsMesh<T_Mesh> &mesh)
  {
    mesh.getPointData().insertNextScalarData(color, "color");
    // adjustment for new ViennaLS
  }

  lsSmartPointer<treeNode> getRoot()
  {
    return root;
  }

  nodes_vector &getTreeNodes() { return treeNodes; }

  /// public interface for the neighborhood query (for a single point)
  template <class Point, class iteratorType = iterator>
  iteratorType getNearestNeighbors(const Point &point)
  {
    // if (not frozen)
    // {
    //   lsMessage::getInstance()
    //       .addWarning("lsTree was not built. Returning empty iterator.")
    //       .print();
    //   return nullptr; //iteratorType(*this);
    // }
    lsSmartPointer<treeNode> foundBin = getBin(point);
    return data[foundBin->leafNum].begin();
  }

  /// Get the nearest Neighbors in form of one or more bins.
  ///
  /// Returns an iterator to the nearest points found.
  /// Since the lsTree partioned the space into multiple bins, the resulting
  /// iterator of this query might refer to multiple disconnected segments
  /// in memory - even if the original container was contingent
  /// (e.g. std::vector, std::array).
  template <class iteratorType = iterator>
  iteratorType getNearestNeighbors(Data_Container &points) const
  {
    if (treeNodes.empty())
    {
      lsMessage::getInstance()
          .addWarning("lsTree was not built. Returning empty iterator.")
          .print();
      return iteratorType(*this);
    }

    nodes_vector found_nodes;
    for (auto &pt : points)
      found_nodes.push_back(getBin(pt));

    return iteratorType(found_nodes, *this);
  }

  /*
   * ---- Configure TREE ---- *
   */

  void setMaxDepth(std::size_t newDepth) { this->maxDepth = newDepth; }

  void setMaxNumBins(std::size_t newLimit) { this->maxNumBins = newLimit; }

  void setMaxPointsPerBin(std::size_t newLimit) { this->maxPointsPerBin = newLimit; }

  void setBestSplitMode(bool newBestSplitMode) { this->useBestSplit = newBestSplitMode; }

  /// Reconfigure the trees parameters.
  ///
  /// Already built lsTree will not be rebuilt unless 'apply()' is called.
  void configureTreeSettings(std::size_t maxDepth_, std::size_t maxNumBins_,
                             std::size_t maxPointsPerBin_)
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
#ifndef NDEBUG // if in debug build
    {
      auto numLeafs = this->endLeafs() - this->beginLeafs();
      bool check = numBins == numLeafs;
      if (!check)
        lsMessage::getInstance().addDebug("lsTree: numBins: " + std::to_string(numBins) + " <= numLeafs:" + std::to_string(numLeafs)).print();
    }
#endif
    return numBins;
  }

  size_type getNumberOfLeafs() const
  {
    return this->getNumberOfBins();
  }

  size_type getNumberOfNodes() const { return treeNodes.size(); }

  uint getDepth() const { return depth; }

  const std::string getTreeType() { return tree_type; }

  /// prints parameters of the tree
  void printInfo()
  {
    std::cout << getTreeType() << std::endl;
    std::cout << "lsTree Settings" << std::endl;
    std::cout << "depth: " << depth << std::endl;
    std::cout << "maxDepth: " << maxDepth << " --> depth: " << depth
              << std::endl;
    std::cout << "maxNumBins: " << maxNumBins << " --> numBins: " << numBins
              << std::endl;
    std::cout << "maxPointsPerBin: " << maxPointsPerBin
              << " --> largestBinSize: " << largestBinSize << std::endl;
  }

  /// prints the tree to console
  void printBT(bool printPoints = false)
  {
    if (frozen)
      printBT("", treeNodes[0], false, printPoints);
  }

  /// dumps the points in the leaf nodes into a python file for easy import and visualization (it's a pretty cheap and dirty way!)
  void dumpBins(std::ostream &os = std::cout) const
  {
    unsigned idx = 0;
    os << "bins = []\n";
    for (auto bin = treeNodes.begin() + startLeafs; bin != treeNodes.end(); ++bin)
    {
      os << "# bin#: " << idx++ << "\n";
      os << "bins.append([";
      printPointsInNode(*bin, os);
      os << "])\n";
    }
    os << std::endl;
  }

  /// prints points in the node to console
  void printPointsInNode(lsSmartPointer<treeNode> node, std::ostream &os = std::cout) const
  {
    auto &bin = data[node->leafNum];
    for (auto &pt : bin)
    {
      lsInternal::print_point(pt, os);
      os << ",";
    }
    os << std::endl;
  }

#pragma Private Methods
private:
  /// helper function that prints a tree node
  void printBT(const std::string &prefix, lsSmartPointer<treeNode> node,
               bool isLeft, bool printPoints)
  {
    if (node != nullptr)
    {
      std::cout << prefix;
      std::cout << (isLeft ? "├──" : "└──");

      if (node->identifier == "")
      {
        std::string dim = "z";
        if (node->dimensionToSplit() == Dim::y)
          dim = "y";
        if (node->dimensionToSplit() == Dim::x)
          dim = "x";
        std::cout << "(" << node->color << " / " << dim << " / " << node->median
                  << ")" << std::endl;
      }
      else
      {
        std::cout << "(" << node->color << " / " << node->identifier << ")"
                  << std::endl;
      }
      if (node->isLeaf && printPoints)
      {
        printPointsInNode(node);
      }

      // enter the next tree level - left and right branch
      printBT(prefix + (isLeft ? "│   " : "    "), node->left, true, printPoints);
      printBT(prefix + (isLeft ? "│   " : "    "), node->right, false, printPoints);
    }
  }

  /*
   * ---- UTILITY METHODS ---- *
   */

  /// it computes the dimension of the best split for a given data collection (given via iterators), i.e. looks for the dimension with the largest span.
  template <typename It>
  unsigned computeNextBestSplitDimension(It begin, It end) const
  {

    std::array<std::array<value_type, 2>, D> limits;
    // fill with (lowest value, highest value) - pairs --> depends on D!
    limits.fill(std::array<value_type, 2>{std::numeric_limits<value_type>::max(), std::numeric_limits<value_type>::lowest()});

    for (auto &idx = begin; idx != end; ++idx)
    {
      auto &point = data[0][*idx];
      for (unsigned dim = 0; dim < D; ++dim)
      {
        limits[dim][0] = (point[dim] < limits[dim][0]) ?: limits[dim][0];
        limits[dim][1] = (point[dim] < limits[dim][1]) ?: limits[dim][1];
      }
    }

    std::array<value_type, D> spans;
    // fill with (lowest value, highest value) - pairs --> depends on D!
    spans.fill(std::numeric_limits<value_type>::min());
    for (unsigned dim = 0; dim < D; ++dim)
    {
      spans[dim] = limits[dim][1] - limits[dim][0];
    }

    auto cmp = [](const value_type &a, const value_type &b)
    {
      return (std::abs(a) < std::abs(b));
    };

    auto result = std::max_element(spans.begin(), spans.end(), cmp);
    unsigned highest = std::distance(spans.begin(), result);

    return highest;
  }

  /// applies the "color" (value between [-numBins, +numBins]) to the nodes.
  /// This information is for Debugging only and has no other use.
  void applyColorToNodes()
  {
    int b0 = numBins;
    // auto &root = treeNodes.front();
    colorChild(root->left, root->color, b0, true);
    colorChild(root->right, root->color, b0, false);
  }

  void colorChild(lsSmartPointer<treeNode> node, size_type parentColor, int b0,
                  bool isLeft)
  {
    auto b = lsInternal::divByPow2(b0, node->level);
    node->color = isLeft ? parentColor - b : parentColor + b;
    if (node->left)
      colorChild(node->left, node->color, b0, true);
    if (node->right)
      colorChild(node->right, node->color, b0, false);
  };

  /// sorts a range of data in a certain dimension
  /// begin, end are iterators over the array of indices for the final ordering!
  template <class It>
  void sortByDim(It begin, It end, uint dim)
  {
    std::sort(begin, end, [this, &dim](const auto &left, const auto &right)
              {
                // return data[0][left][dim] > data[0][right][dim];
                return data[0][left][dim] < data[0][right][dim];
              });
  }

  /// Sorts a vector of 3D points (point cloud) by generating a sorted
  /// vector of indices (of length size) and move-returns it.
  /// Does NOT actually sort the vector (in place)!
  template <class VectorIt>
  std::vector<size_type> argsortInDimension(VectorIt toSortBegin, size_type dim,
                                            size_type size)
  {
    std::vector<size_type> result(size);
    std::generate(result.begin(), result.end(),
                  [n = 0]() mutable
                  { return n++; });
    std::sort(
        result.begin(), result.end(),
        [toSortBegin,
         dim](auto left,
              auto right) { // left, right : size_type -> indices of vector toSort
          return toSortBegin[left][dim] <
                 toSortBegin[right][dim]; // sort in ascending order
        });
    return std::move(result);
  }

  /// Sorts a vector data in place by ordering given via vector indices (same as
  /// sortByIndices but in place)
  template <class Vector, class VectorIndices>
  void sortByIndicesInplace(Vector &data, VectorIndices &indices)
  {
    auto beginOfLeafs = this->beginLeafs();
    auto endOfLeafs = this->endLeafs();
    auto numberOfLeafs = endOfLeafs - beginOfLeafs; // == number of bins
    auto data0 = data[0];
#ifndef NDEBUG // if in debug build
    {
      bool check = data0.size() == indices.size();
      if (!check)
        lsMessage::getInstance().addDebug("lsTree: data.size(): " + std::to_string(data0.size()) + " <= indices.size():" + std::to_string(indices.size())).print();

      check = numberOfLeafs == numBins;
      if (!check)
        lsMessage::getInstance()
            .addDebug("lsTree: numberOfLeafs: " + std::to_string(numberOfLeafs) + " =? " + std::to_string(numBins))
            .print();
    }
#endif
    // Allocating the vector we will copy into
    //                    #bins          each bin has #lBS number of points (=std::array<double, 3>)
    Data_Container result(numberOfLeafs, Bin(largestBinSize, point_type{0, 0, 0}));
    auto binIdx = 0;
    for (auto nodeIt = beginOfLeafs; nodeIt != endOfLeafs; ++nodeIt, ++binIdx)
    {
      lsSmartPointer<treeNode> current_node = *nodeIt;
      auto start = current_node->start;
      auto stop = current_node->stop; // get number of points in node

#ifndef NDEBUG // if in debug build
      {
        bool check = binIdx <= numBins;
        if (!check)
          lsMessage::getInstance().addDebug("lsTree: binIdx: " + std::to_string(binIdx) + " <= " + std::to_string(numBins) + " --> [NO]").print();

        const auto size = stop - start;
        check = size <= largestBinSize;
        if (!check)
          lsMessage::getInstance().addDebug("lsTree: newIdx: " + std::to_string(size) + " <= " + std::to_string(largestBinSize) + " --> [NO]").print();
      }
#endif
      current_node->leafNum = binIdx;

      for (unsigned i = start, newIdx = 0; i < stop; ++i, ++newIdx)
      {
        auto indexToMove = indices[i];
        auto ptToMove = data0[indexToMove];
        result[binIdx][newIdx] = ptToMove;
      }
    }
    data = std::move(result);
  }
}; // end class lsTree

#endif // LS_TREE
