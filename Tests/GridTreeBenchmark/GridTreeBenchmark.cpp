// STL
// #include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
// #include <random>
#include <vector>
#include <chrono>
#include <utility>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToDiskMesh.hpp>
//#define LS_TREE
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>
#include <lsTestAs>

using results_container = std::vector<std::tuple<int, int, double>>;

void print_results(results_container results)
{
  std::cout << "D   |   N   |  t[ms]\n";
  std::cout << "--------------------\n";
  for (const auto &row : results)
  {
    auto [D, N, t] = row;
    // auto D = std::get<0>(row);
    // auto N = std::get<1>(row);
    // auto t = std::get<2>(row);
    std::cout << std::setw(4) << D << "|" << std::setw(7) << N << "|" << std::setw(7) << t << std::endl;
  }
};

template <int D = 2>
const int benchmark()
{
  int before = omp_get_max_threads();
  omp_set_num_threads(4);
  std::cout << "Using " << omp_get_max_threads() << " (of max " << before << ") threads for this test." << std::endl;

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  constexpr double MAX_RADIUS = 6.1, MIN_RADIUS = 1;
  results_container results; //int(MAX_RADIUS - MIN_RADIUS + 1));

  for (double radius = 1; radius < MAX_RADIUS; ++radius)
  {
    // radius = 1 --> 25 po{ints
    // radius = 2 --> 86 points‚
    // radius = 3 --> 230 points
    // radius = 4 -->  points
    const hrleVectorType<double, D> centre(5., 0., 0.);

    lsMakeGeometry<double, D>(
        levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
        .apply();
    lsToMesh<double, D>(levelSet, mesh).apply();

    auto N = levelSet->getDomain().getNumberOfPoints();

    // const std::string file_prefix = "/Users/peterholzner/Code/ViennaTools/ViennaLS/lsTree/Tests/GridTreeBenchmark/meshes/" + std::to_string(D) + "D_N" + std::to_string(N) + "_";
    const std::string file_prefix = "./meshes/" + std::to_string(D) + "D_N" + std::to_string(N) + "_";
    lsVTKWriter<double>(mesh, file_prefix + "Mesh.vtk").apply();

    // ---- Benchmark
    std::cout << "D: " << D;
    std::cout << "\nNumber of points: " << N << std::endl;
    std::cout << "Building tree... " << std::flush;
    const auto start = std::chrono::high_resolution_clock::now();
    auto tree = lsTree<double, D>(mesh);
    tree.apply();
    const auto stop = std::chrono::high_resolution_clock::now();
    const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    results.push_back(std::make_tuple(D, N, time));
    std::cout << "DONE!";
    std::cout << "\nt: " << time << "ms";
    std::cout << "\n---------------------------------------------------" << std::endl;

    lsVTKWriter(mesh, lsFileFormatEnum::VTU, file_prefix + "TreeMesh.vtu")
        .apply();
    tree.printInfo();
    // tree.printTree();
    // tree.printTreeByLevel();
    // tree.printBT();
  }
  std::cout << "\n---------------------------------------------------" << std::endl;
  print_results(results);
  return EXIT_SUCCESS;
};

int main(int argc, char **argv)
{
  int success_2d_test = EXIT_SUCCESS, success_3d_test = EXIT_SUCCESS;
  std::cout << "\n###################################################";
  std::cout << "\n######################## D = 2 ####################";
  std::cout << "\n###################################################" << std::endl;
  success_2d_test = benchmark<2>();
  std::cout << "\n###################################################";
  std::cout << "\n######################## D = 3 ####################";
  std::cout << "\n###################################################" << std::endl;
  success_3d_test = benchmark<3>();
  if (success_2d_test == EXIT_SUCCESS && success_3d_test == EXIT_SUCCESS)
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return 0;
}
