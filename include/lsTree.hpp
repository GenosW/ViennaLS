//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <numeric>
//#include<cstdint>
//#include <iterator>
#include <lsDimension.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>
//#include <utility>
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
  constexpr int64_t pow2(int exp) { return 1 << exp; }

  constexpr int64_t divByPow2(int num, int exp) { return num >> exp; }

  template <int D>
  std::size_t inline convertToIndex(lsInternal::Dimension<D> dim)
  {
    return dim.toArrayIndex();
  }

  enum struct meanEnum : unsigned
  {
    MEDIAN = 0,
    ARITHMETIC = 1,
    // GEOMETRIC = 2,
    HARMONIC = 3,
  };

  template <typename Point>
  void print_point(Point &p)
  {
    std::cout << "[" << p[0] << ", " << p[1] << ", " << p[2] << "]";
    // std::cout << "(" << *p[0][0] << "," << *p[0][1] << "," << *p[0][2] << ")" << std::endl;
  };
} // end namespace lsInternal

/// Tree structure that seperates a domain
///
/// TODO: clean up parameters, member functions and attributes
/// TODO: order
///
/// Types:
///
/// using value_type = T;
/// using size_type = std::size_type
/// using point_type = std::array<T, 3>;
/// using Data_Container = std::vector<point_type>;
/// using iterator = treeIterator<point_type>;
/// using const_iterator = constTreeIterator<point_type>;
/// using range_vector = typename iterator::range_vector;
/// using nodes_vector = std::vector<lsSmartPointer<treeNode>>;
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
      // os << node.color << ": [" << node.start << ", " << node.top << "]";
      os << "Hi";
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
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  size_type maxDepth = 4; // TODO: maxDepth prop log10(N)
  size_type maxNumBins = lsInternal::pow2(maxDepth);
  size_type maxPointsPerBin = 4;
  lsInternal::meanEnum meanMethod = lsInternal::meanEnum::MEDIAN;

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
  nodes_iterator beginLeafs(void)
  {
    return std::next(treeNodes.begin(), startLeafs);
  }

  nodes_const_iterator beginLeafs(void) const
  {
    return std::next(treeNodes.begin(), startLeafs);
  }

  // typename nodes_vector::iterator
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

  //
  // lsTree(ContainerToWrap container)
  //     : data(container.getData()),
  //       N(container.getData().size()) {}

  //~lsTree() {}

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
      root = lsSmartPointer<treeNode>::New(0, D, 0, N, data[0][sortedPoints[N / 2]][D - 1]);

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
  /// calculates index of level based dimensionality of lsTree instance.
  /// [DEPRECATED]
  size_type constexpr getNextDimIdx(size_type level)
  {
    return D - 1 - level % (D);
  }

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

      // Sort Range of the root
      // if (level > 1) // otherwise already sorted
      Dim nextSplitDim = root->dimensionToSplit--;
      // auto sortDimIdx = lsInternal::convertToIndex(root->dimensionToSplit);
      auto sortDimIdx = lsInternal::convertToIndex(nextSplitDim);
      // sortByDim(sortedPoints.begin() + start, sortedPoints.begin() + stop,
      //           sortDimIdx);

      // Make 2 new nodes
      auto leftMedianIndex = start + leftRange / 2;
      auto rightMedianIndex = start + leftRange + rightRange / 2;

      T leftAverage = 0, rightAverage = 0;
      if (meanMethod == lsInternal::meanEnum::MEDIAN)
      {
        // median = lookup Point from original data
        auto get_median = [this, &sortedPoints, sortDimIdx](size_type medianIndex, size_type range)
        {
          auto result = (range % 2 == 0) ? (
                                               (data[0][sortedPoints[medianIndex]][sortDimIdx] +
                                                data[0][sortedPoints[medianIndex - 1]][sortDimIdx]) /
                                               2)
                                         : (data[0][sortedPoints[medianIndex]][sortDimIdx]);
          return result;
        };
        // T median = data[0][sortedPoints[leftMedianIndex]][sortDimIdx];
        // if (leftRange % 2 == 0)
        // {
        //   median = (median + data[0][sortedPoints[leftMedianIndex - 1]][sortDimIdx]) / 2;
        // }
        // arithmetic middle if an even number of elements are in the leftRange
        // arithmetic middle if an even number of elements are in the rightRange
        // median = (rightRange % 2 == 0) ? ((data[0][sortedPoints[rightMedianIndex]][sortDimIdx] + data[0][sortedPoints[rightMedianIndex - 1]][sortDimIdx]) / 2) : data[0][sortedPoints[rightMedianIndex]][sortDimIdx];
        leftAverage = get_median(leftMedianIndex, leftRange);
        rightAverage = get_median(rightMedianIndex, rightRange);
      }
      else if (meanMethod == lsInternal::meanEnum::ARITHMETIC)
      {
        // std::reduce(std::execution::par, v.cbegin(), v.cend())

        for (uint idx = start; idx < start + leftRange; ++idx)
        {
          leftAverage += data[0][sortedPoints[idx]][sortDimIdx];
        }
        leftAverage /= leftRange;

        for (uint idx = start + leftRange; idx < stop; ++idx)
        {
          rightAverage += data[0][sortedPoints[idx]][sortDimIdx];
        }
        rightAverage /= rightRange;
      }
      else if (meanMethod == lsInternal::meanEnum::HARMONIC)
      {
        for (uint idx = start; idx < start + leftRange; ++idx)
        {
          leftAverage += 1 / data[0][sortedPoints[idx]][sortDimIdx];
        }
        leftAverage = leftRange / leftAverage;

        for (uint idx = start + leftRange; idx < stop; ++idx)
        {
          rightAverage += 1 / data[0][sortedPoints[idx]][sortDimIdx];
        }
        rightAverage = rightRange / rightAverage;
      }
      else
      {
        leftAverage = 0;
        rightAverage = 0;
      }

      root->left = lsSmartPointer<treeNode>::New(level, nextSplitDim, start, start + leftRange, leftAverage, root);
      root->right = lsSmartPointer<treeNode>::New(
          level, nextSplitDim, start + leftRange, stop, rightAverage, root);

      sortByDim(sortedPoints.begin() + start, sortedPoints.begin() + stop,
                sortDimIdx);

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

  /// builds a treeNode of the tree and returns a pointer
  ///
  /// checks if it should be a leaf treeNode
  /// TODO: Messy template
  // template <class Vector, class VectorData>
  // lsSmartPointer<treeNode> buildByDepth(Vector &data, VectorData &originalData,
  //                                       size_type start, size_type stop,
  //                                       size_type level)
  // {
  //   size_type range = stop - start;
  //   if ((range < maxPointsPerBin) || (level > maxDepth))
  //     return nullptr;
  //   size_type leftRange = range / 2;
  //   depth = (depth < level) ? level : depth;
  //   // size_type thisIndex = index++;

  //   // thisNodes split dimension
  //   size_type dim = D - level % (D + 1); // offset
  //   //++index;

  //   T median = originalData[data[0][start + leftRange / 2]][dim];
  //   auto thisNode =
  //       lsSmartPointer<treeNode>::New(level, dim, start, stop, median);
  //   treeNodes.push_back(thisNode);
  //   // thisNode->setRange(start, stop);
  //   // thisNode->level = level;
  //   // thisNode->dimensionToSplit = dim;

  //   thisNode->left =
  //       buildByDepth(data, originalData, start, start + leftRange, level + 1);
  //   thisNode->right =
  //       buildByDepth(data, originalData, start + leftRange, stop, level + 1);
  //   // thisNode->isLeaf = (thisNode->left != nullptr || thisNode->right !=
  //   // nullptr ) ? false : true;

  //   return thisNode;
  // }

  //---- QUERY METHODS ----

  // private: // TODO: make it private again
public:
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

  // iteratorType getNearestNeighbors(const point_type &point)
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

  void setMeanMethod(lsInternal::meanEnum newMethod)
  {
    meanMethod = newMethod;
  }

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

  /// prints parameters
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

  // void printBT(const std::string &prefix, treeNode &node,
  //              bool isLeft, printPoints)
  // {
  //   if (frozen)
  //   {
  //     lsSmartPointer<treeNode> nodePointer(&node);
  //     printBT(nodePointer)
  //     // std::cout << prefix;
  //     // std::cout << (isLeft ? "├──" : "└──");

  //     // if (node.identifier == "")
  //     // {
  //     //   std::string dim = "z";
  //     //   if (node.dimensionToSplit() == Dim::y)
  //     //     dim = "y";
  //     //   if (node.dimensionToSplit() == Dim::x)
  //     //     dim = "x";
  //     //   std::cout << "(" << node.color << " / " << dim << " / " << node.median
  //     //             << ")" << std::endl;
  //     // }
  //     // else
  //     // {
  //     //   std::cout << "(" << node.color << " / " << node.identifier << ")"
  //     //             << std::endl;
  //     // }

  //     // // enter the next tree level - left and right branch
  //     // printBT(prefix + (isLeft ? "│   " : "    "), node.left, true);
  //     // printBT(prefix + (isLeft ? "│   " : "    "), node.right, false);
  //   }
  // }

  void printBT(bool printPoints = false)
  {
    if (frozen)
      printBT("", treeNodes[0], false, printPoints);
  }

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

  void printPointsInNode(lsSmartPointer<treeNode> node, std::ostream &os = std::cout) const
  {
    auto &bin = data[node->leafNum];
    for (auto &pt : bin)
    {
      lsInternal::print_point(pt);
      os << ",";
    }
    os << std::endl;
  }

#pragma Private Methods
private:
  /*
   * ---- UTILITY METHODS ---- *
   */

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
    std::sort(begin, end, [beginSortBy, endSortBy](auto left, auto right)
              {
                // sort in ascending order
                for (auto iterator_to_sorted_data = beginSortBy;
                     iterator_to_sorted_data < endSortBy; ++iterator_to_sorted_data)
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
    std::sort(begin, end, [this, &dim](const auto &left, const auto &right)
              {
                return data[0][left][dim] < data[0][right][dim];
                // false...switch order of left and right });
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

private:
  inline size_type getRoot(size_type treeNode)
  {
    return std::ceil((double)treeNode / 2) - (1);
  }

  inline size_type getLeftChild(size_type treeNode) { return treeNode * 2 + 1; }

  inline size_type getRightChild(size_type treeNode) { return treeNode * 2 + 2; }
}; // end class lsTree

#endif // LS_TREE
