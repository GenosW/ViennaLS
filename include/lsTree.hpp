//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <algorithm>
#include <iterator>
#include <vector>
#include <lsToDiskMesh.hpp>
#include <lsSmartPointer.hpp>

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
///
///
/// Maximum depth is 4?
template <class T, int D>
class lsTree
{
public:
  using point_type = std::array<T, 3>;
  typedef typename std::array<T, 3> point_type2;

private:
  // The mesh we're building a tree for
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  lsSmartPointer<lsMesh<T>> newMesh = nullptr;

  // PARAMETERS
  std::string tree_type = "kd-Tree";
  int depth = 0;
  int numBins = 1; // numBins = 2^depth
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  int maxDepth = 4;
  int maxNumBins = ipower(2, maxDepth);
  int maxPointsPerBin = 10;

  struct node
  {
    // std::vector<point_type> points;
    size_t start = 0;
    size_t stop = 0;
    size_t level = 0;
    size_t size = 0;
    point_type center;
    point_type minimumExtent;
    point_type maximumExtent;
    lsSmartPointer<node> left = nullptr;
    lsSmartPointer<node> right = nullptr;

    // node(const point_type& ct, const point_type& tl, const point_type& br) : center(ct), extent({tl, br}) left(nullptr), right(nullptr) {
    // }
    node()
    {
    }

    node(const point_type &ct) : center(ct)
    {
    }

    node(const point_type &ct, const point_type &tl, const point_type &br) : center(ct)
    {
    }

    void setRange(size_t newStart, size_t newStop)
    {
      start = newStart;
      stop = newStop;
      size = stop - start;
    }

    size_t getSize()
    {
      return size;
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

    void split()
    {
      left = lsSmartPointer<node>::New();
      size_t med = (stop - start) / 2;
      left->setRange(start, med);
      right = lsSmartPointer<node>::New();
      right->setRange(med, stop);
    }
  };
  lsSmartPointer<node> root = nullptr;
  std::vector<lsSmartPointer<node>> nodes;

  // Methods
  // void build() {}

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

    std::cout << "Partitioning... " << std::endl;
    // std::vector<point_type> below(begin, begin + medianPos);
    // std::vector<point_type> above(begin + medianPos, end);
    std::vector<T> below(medianPos, 0.);
    std::vector<T> above(medianPos, 1.);
    std::cout << "below.size: " << below.size() << std::endl;
    std::cout << "above.size: " << above.size() << std::endl;
    // mesh->insertNextVectorData(below, "below");
    // mesh->insertNextVectorData(above, "above");
    mesh->insertNextScalarData(below, "below");
    mesh->insertNextScalarData(above, "above");

    size_t index = 0;
    root = build(begin, 0, N, 0, index);
  }

  template <class VectorIt>
  lsSmartPointer<node> build(VectorIt begin, size_t start, size_t stop, size_t level, size_t &index)
  {
    size_t size = stop - start;
    if (size < maxPointsPerBin)
      return nullptr;
    if (level > maxDepth)
      return nullptr;
    size_t halfSize = size / 2;
    // std::nth_element(begin, begin+halfSize, begin+size); // doesn't work cause we have triples here!
    // partitionInDimension(dim, begin+start, begin+stop, begin[halfSize]);
    size_t nextLevel = (level + 1);
    depth = (depth < level) ? level : depth;
    size_t thisIndex = index;
    index += 1;
    auto thisNode = lsSmartPointer<node>::New();
    nodes.push_back(thisNode);
    thisNode->setRange(start, stop);
    thisNode->level = level;
    thisNode->left = build(begin, start, start + halfSize, nextLevel, index);
    thisNode->right = build(begin, start + halfSize, stop, nextLevel, index);

    return thisNode;
  }

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

  // template <class VectorIt>
  // size_t partitionInDimension(size_t dim, VectorIt begin, VectorIt end, size_t median)
  // {
  //   std::partition(begin, end, [dim, median](const auto &pos) { return pos[dim] < median; })
  //   return std::distance(begin, end) / 2
  // }

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

  void printInfo()
  {
    std::cout << getTreeType() << std::endl;
    std::cout << "depth: " << depth << std::endl;
    std::cout << "numBins: " << numBins << std::endl;
    std::cout << "maxDepth: " << maxDepth << std::endl;
    std::cout << "maxNumBins: " << maxNumBins << std::endl;
    std::cout << "maxPointsPerBin: " << maxPointsPerBin << std::endl;
  };

  void printTree()
  {
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      std::string leaf = (nodes[i]->left == nullptr && nodes[i]->right == nullptr) ? "#" : "";
      std::cout << "node " << std::setw(3) << i << "(L" << nodes[i]->level << "): [" << nodes[i]->start << ", " << nodes[i]->stop << ")" << leaf << std::endl;
    }
  };

  void printTree2()
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
