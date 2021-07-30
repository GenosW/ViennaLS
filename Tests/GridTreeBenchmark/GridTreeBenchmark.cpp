// STL
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <utility>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
// #include <lsToDiskMesh.hpp>
#include <lsTest.hpp>
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>
#include <lsTestAsserts.hpp>

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

void save_results(results_container results, const std::string path)
{
  std::ofstream file;
  file.open(path);
  file << "D,N,t\n";
  for (const auto &row : results)
  {
    auto [D, N, t] = row;
    // auto D = std::get<0>(row);
    // auto N = std::get<1>(row);
    // auto t = std::get<2>(row);
    file << D << "," << N << "," << t << std::endl;
  }
  file.close();
};

template <int D>
const int benchmark()
{
  setup_omp(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  constexpr double MAX_RADIUS = 20.1, MIN_RADIUS = 1;
  results_container results; //int(MAX_RADIUS - MIN_RADIUS + 1));

  for (double radius = 1; radius < MAX_RADIUS; radius += 2)
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
    // Get reference geometry parameters
    auto N = levelSet->getDomain().getNumberOfPoints();
    LSTEST_ASSERT(mesh->getNodes().size() == N);

    const std::string file_prefix = "./meshes/" + std::to_string(D) + "D_N" + std::to_string(N) + "_";
    lsVTKWriter<double>(mesh, file_prefix + "Mesh.vtk").apply();

    // ---- Benchmark
    std::cout << "D: " << D;
    std::cout << "\nNumber of points: " << N << std::endl;
    std::cout << "Building tree... " << std::flush;
    const auto start = std::chrono::high_resolution_clock::now();
    auto tree = lsTree<double, D>(mesh).apply();
    const auto stop = std::chrono::high_resolution_clock::now();
    const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    results.push_back(std::make_tuple(D, N, time));
    std::cout << "DONE!";
    std::cout << "\nt: " << time << "ms";
    std::cout << "\n---------------------------------------------------" << std::endl;

    lsVTKWriter(mesh, lsFileFormatEnum::VTU, file_prefix + "TreeMesh.vtu")
        .apply();
    // tree.printInfo();
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
  lsTest test_2d = INIT_LSTEST(test_2d);
  MAKE_LSTEST(test_3d);

  std::cout << "\n###################################################"
            << "\n######################## D = 2 ####################"
            << "\n###################################################"
            << std::endl;
  test_2d.run(benchmark<2>);
  std::cout << "\n###################################################"
            << "\n######################## D = 3 ####################"
            << "\n###################################################"
            << std::endl;
  test_3d.run(benchmark<3>);

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  check(test_2d, test_3d);

  if (all(test_2d.wasSuccess(), test_3d.wasSuccess()))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return 0;
}
