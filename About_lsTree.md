# lsTree

The goal of this fork is to introduce a tree data structure that provides a nearest neighbor query interface. The lsTree data structure should build a BST (binary search tree), reorder or copy the input data container in such a way that the result of the neighbor query is:

- one or multiple iterators over a collection of points ("bin" or "bucket") that represent a reasonable neighborhood of the queried point,
- each bin should be laid out in a way to provide optimal access and iteration times, i.e. it should be contingent in memory in order to be fed to the parallelized algorithms in this library.

Although other tree structure were considered, and are likely still viable such as Quad-/Octrees, we chose to test the implementation of a kd-Tree type structure with bins (collections of points) in the leaf nodes.

## kdTree

- [kd-Tree on Wikipedia](https://en.wikipedia.org/wiki/K-d_tree)
- [Rosetta Code kd-Tree implementation](https://rosettacode.org/wiki/K-d_tree)

A kd-Tree adaptively divides a collection of data (points). For non-leaf nodes, the split is based on the median value of the respective collection of points that the node represents. In the case of multidimensional data, such as is usually the case here, one of the dimensions has to be chosen to compute the median (each splitting plane is parallel to one of the axis). The choice can be either (1) at random, (2) follow a set cycle (e.g. z->y->x->z->...) or (3) can be determined via some other measure such as choosing the dimension with the biggest span.

(2) and (3) are already implemented. (2) will also lead each node of the same level/depth to have the same splitting dimension, where as (1) and (3) might result in different splitting dimensions for each node (in a level).

## Parameters

The lsTree data structure can be adjusted via certain parameters settings that control various features of the tree, such as the maximum depth (`maxDepth`) or maximum number of points per bin (`maxPointsPerBin`). These can be set via setter methods.

Settings:

- maxDepth = 4
- maxNumBins = lsInternal::pow2(maxDepth)
- maxPointsPerBin = 4

The actual features of the tree and then recorded in these attributes:

- depth = 0;
- N = 0;       // number of points in data / size of data
- numBins = 1; // numBins = 2^depth
- largestBinSize = 0;

## Conclusion so far

The kd-tree implementation does not work yet, although the reasons are not clear to me yet. For some test geometries it behaves as expected for others it doesn't.

The unit tests to verify this are located in the folder `GridTree`. Other tests associated with this fork are in the folders `GridTreeIterator`, `GridTreeNeighbors` and `Dimensions`, although some of them correspond to earlier versions and are now defunct!

Additional notes, that were made during the implementation, can be found in the folder `lsTree_notes`.
