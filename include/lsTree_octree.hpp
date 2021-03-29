//#define LS_TREE
#ifndef LS_TREE
#define LS_TREE
#include <vector>
#include <lsToDiskMesh.hpp>
#include <math.h>

// (x > 0) : (x < 0) -> Close : Far
// (y > 0) : (y < 0) -> Right : Left
// (z > 0) : (z < 0) -> Top : Bottom
// The direction/Reihenfolge of these numbers might need tweaking!
#define CRT 0    // close right top <--> (x>0) && (y>0) && (z>0)
#define FRT 1    // far right top   <--> (x<0) && (y>0) && (z>0)
#define FLT 2    // far left top    <--> (x<0) && (y<0) && (z>0)
#define CLT 3    // close left top  <--> (x>0) && (y<0) && (z>0)
#define CRB 4    // close right bot <--> (x>0) && (y>0) && (z<0)
#define FRB 5    // far right bot   <--> (x<0) && (y>0) && (z<0)
#define FLB 6    // far left bot    <--> (x<0) && (y<0) && (z<0)
#define CLB 7    // close left bot  <--> (x>0) && (y<0) && (z<0)

inline int ipower(int N, int exp) {return (exp>1) ? N*ipower(N, exp-1) : N;};

// template <int Dim=3> struct domain {
//   std::vector<int> data;

//   domain(int N){
//     size_t S = ipower(N, Dim);
//     data.reserve(S);
//     int start = 0;
//     auto next_value = [&start] () {return start++;};
//     std::generate(data.begin(), data.end(), next_value);
//   }

//   size_t size() {
//     return data.size();
//   }
// };

/// Tree structure that seperates a domain
/// Tree type is based on dimension:
///   - D=2... Quadtree
///   - D=3... Octree
/// Maximum depth is 4?
template <class T, int D> class lsTree { 
  lsSmartPointer<lsToDiskMesh> mesh = nullptr;
  lsSmartPointer<lsDomain<T, D>> domain = nullptr;

  // PARAMETERS
  // TODO: Validate the maxDepth setting. Too deep/too shallow/just right?
  int maxDepth = 4;
  // TODO: Rectangular grid???
  
  int pointsPerDim = 20;
  // TODO: (BR) Reduce number of parameters to minimum needed.
  // They are here to have something to look at and have default values.
  // Remove before 'release'
  int binSize = ipower(20, D);
  int numChildren = ipower(2, D);
  int minSizePerDirection = pointsPerDim * 2;
  int minSizeGrid = numChildren * binSize;

  // lsTree node
  struct Node {
    // 
    //coords boundary;
    // Children
    // Vector of trees for scalability
    // Access via integer index or DEFINES above
    std::vector<lsSmartPointer<lsPointData>> points; // length = numChildren
    std::vector<lsSmartPointer<Node>> children; // length = numChildren
    bool hasChildren;

    void split() {
      hasChildren = true:
      children.assign(numChildren);
    }

    void makeLeaf() {
      hasChildren = false;
    }

    void storePoint() {
      // Add to points
    }

    bool contains(lsPointData point) {
      // Check if Node contains the point based on boundary
      return false;
    }
  }
  Node* root; // TODO: (BR) Convert to proper pointer

  // Methods
  // void build() {}
  // bool split() {}
  // bool insert(point) {}
  // region query(point) {}
  // regions query(boundary) {}

public:
  lsTree() {}

  lsTree(lsSmartPointer<lsDomain<T, D>> passedDomain,
               lsSmartPointer<lsToDiskMesh> passedMesh)
      : domain(passedDomain), mesh(passedMesh) {
      }

  //~lsTree() {}

  /// static build --> no need to restructure since grid doesn't change?
  void build() {
    
    // TODO: Implement build

    // 1. Check how many points the current node represents
    // 2. If too large, split into quadrants
    // 2. 
  }

  bool insert(lsPointData point) {
    //
  }

  /*
  * ---- DEBUG FUNCTIONS ---- *
  */
  const char* getTreeType() {
    return (D == 2) ? "Quadtree" : "Octree";
  }

  void printInfo(){
      std::cout << getTreeType() << std::endl;
      std::cout << "binSize: " << binSize << std::endl;
      std::cout << "numChildren: " << numChildren << std::endl;
      std::cout << "minSizePerDirection: " << minSizePerDirection << std::endl;
      std::cout << "minSizeGrid: " << minSizeGrid << std::endl << std::endl;
  };

};

#endif // LS_TREE
