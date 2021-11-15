# Notes

## Feedback

- [] unit tests: ctest

## Trees

Interesting stackoverflow about a very similar problem: <https://datascience.stackexchange.com/questions/41345/is-there-any-way-of-ordering-sorting-vectors>

### Quad-/Octree

[C++ octree implementation](https://stackoverflow.com/questions/5963954/fast-templated-c-octree-implementation)

[Octree basics](https://iq.opengenus.org/octree/)

### kd-Tree

## Dev Notes

- Look at Tests as examples more
  - MakeSphere: lsToDiskMesh (modified)
- lsTree will be 'similar' to lsToDiskMesh and lsMesh:
  - std::vector<std::array<T, 3>> nodes: coordinates of (surface) points

## treeIterator

[Iterators](https://www.cplusplus.com/reference/iterator/)

We generally assume that the original data container will be some form of container that is contigous in memory (e.g. std::vector, std::array, etc.) - that just (usually) makes the most sense for simulations.

If you look for just the neighborhood of a singular point, you will simply be given an iterator over a slice of the original data array. However, if you query for a region (e.g. bounding box), you might get more than just "one slice of the original container" - i.e. the resulting iterator might iterate over disconnected segments of the original container.

Maybe the following [stdlib-containers](https://en.cppreference.com/w/cpp/container) and their iterators could provide useful references:

- <https://en.cppreference.com/w/cpp/container/unordered_map>: [begin](https://en.cppreference.com/w/cpp/container/unordered_map/begin) & [end](https://en.cppreference.com/w/cpp/container/unordered_map/end)

- <https://en.cppreference.com/w/cpp/container/unordered_set>: [begin](https://en.cppreference.com/w/cpp/container/unordered_set/begin) & [end](https://en.cppreference.com/w/cpp/container/unordered_set/end)

These containers also operate on data sorted into "buckets" (bins).

*Iterator type*: [std::bidirectional_iterator_tag](https://www.cplusplus.com/reference/iterator/BidirectionalIterator/)

### (Potentially) Non-contiguous sequence of buckets

For now, the treeIterator works with a set (std::vector) of range-pairs:

```cpp
  range_vector ranges;
  // ranges = [(startBin1, stopBin1), (startBin2, stopBin2),
  //           (startBin3, stopBin3), (startBin4, stopBin4), ... ]
```

However, there is a new concept/view-object in the C++20 standard called [std::span](https://en.cppreference.com/w/cpp/container/span) (more info/discussion [here on stackoverflow](https://stackoverflow.com/questions/45723819/what-is-a-span-and-when-should-i-use-one)).
Using std::span does have some practical problems though, mostly to do with the fact that it either requires explicit c++20 support from the compiler or other workaround to make it c++11/14/17 compatible.
That could add some additional complexity to the build process and could make it more difficult for others/users to build (and contribute) OR cause some more trouble to maintain it.
