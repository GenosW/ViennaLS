// STL
#include <chrono>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

// ViennaHRLE
// ViennaLS
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
// #include <lsToDiskMesh.hpp>
#include <lsTest.hpp>
#include <lsTestAsserts.hpp>
#include <lsTree.hpp>
#include <lsVTKWriter.hpp>

#define LINE_SEPARATOR                                                 \
  std::cout << "\n---------------------------------------------------" \
            << std::endl;

using results_container = std::vector<std::tuple<int, double, int, double>>;

void print_results(results_container results, std::string header)
{
  std::cout << "D   |   N   |   r   |  t[ms]";
  std::cout << "\n--------------------\n";
  for (const auto &row : results)
  {
    auto [D, N, r, t] = row;
    // auto D = std::get<0>(row);
    // auto N = std::get<1>(row);
    // auto t = std::get<2>(row);
    std::cout << std::setw(4) << D << "|" << std::setw(7) << N << "|"
              << std::setw(7) << r << "|" << std::setw(7) << t << std::endl;
  }
};

void print_results(results_container results)
{
  print_results(results, "D   |   N   |   r   |  t[ms]");
};

void save_results(results_container results, const std::string path)
{
  std::ofstream file;
  file.open(path);
  file << "D,N,t\n";
  for (const auto &row : results)
  {
    auto [D, N, r, t] = row;
    file << D << "," << N << "," << r << "," << t << std::endl;
  }
  file.close();
};

template <int D>
const int benchmark_cubes()
{
  setup_omp(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  results_container results;

  std::vector<double> testParameters{1, 7, 25, 100, 500, 1000};
  for (auto testParameter : testParameters)
  {
    std::cout << "D: " << D << std::endl;

    // Setup Geometry
    // TODO: make something cubic
    auto length = testParameter;
    std::cout << "length: " << length << std::endl;
    const hrleVectorType<double, D> min(0., 0., 0.), max(length, length, length);
    lsMakeGeometry<double, D>(
        levelSet, lsSmartPointer<lsBox<double, D>>::New(min, max))
        .apply();
    lsToMesh<double, D>(levelSet, mesh).apply();
    //

    // Get reference geometry parameters
    auto N = levelSet->getDomain().getNumberOfPoints();
    LSTEST_ASSERT(mesh->getNodes().size() == N);
    std::cout << "\nNumber of points: " << N << std::endl;

    const std::string file_prefix =
        "./meshes/" + std::to_string(D) + "D_N" + std::to_string(N) + "_";
    lsVTKWriter<double>(mesh, file_prefix + "Mesh.vtk").apply();

    // ---- Benchmark
    std::cout << "Building tree... " << std::flush;
    const auto start = std::chrono::high_resolution_clock::now(); // --- start
    auto tree = lsTree<double, D>(mesh).apply();
    const auto stop = std::chrono::high_resolution_clock::now(); // --- stop

    const auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
            .count();
    results.push_back(std::make_tuple(D, N, radius, time));

    std::cout << "DONE!";
    std::cout << "\nt: " << time << "ms" << std::endl;
    LINE_SEPARATOR;

    lsVTKWriter(mesh, lsFileFormatEnum::VTU, file_prefix + "TreeMesh.vtu")
        .apply();
  }
  LINE_SEPARATOR;
  print_results(results);
  return EXIT_SUCCESS;
};

template <int D>
const int benchmark_sphere()
{
  setup_omp(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  results_container results;

  std::vector<double> testParameters{1, 7, 25, 100, 500, 1000, 2000, 5000};

  for (auto testParameter : testParameters)
  {
    std::cout << "D: " << D << std::endl;

    auto radius = testParameter;
    std::cout << "radius: " << radius << std::endl;
    const hrleVectorType<double, D> centre(5., 0., 0.);
    lsMakeGeometry<double, D>(
        levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
        .apply();
    lsToMesh<double, D>(levelSet, mesh).apply();

    // Get reference geometry parameters
    auto N = levelSet->getDomain().getNumberOfPoints();
    LSTEST_ASSERT(mesh->getNodes().size() == N);
    std::cout << "\nNumber of points: " << N << std::endl;

    const std::string file_prefix =
        "./meshes/" + std::to_string(D) + "D_N" + std::to_string(N) + "_";
    lsVTKWriter<double>(mesh, file_prefix + "Mesh.vtk").apply();

    // ---- Benchmark
    std::cout << "Building tree... " << std::flush;
    const auto start = std::chrono::high_resolution_clock::now(); // --- start
    auto tree = lsTree<double, D>(mesh);
    if (N > 100000)
      tree.setMaxDepth(8);
    tree.apply();
    const auto stop = std::chrono::high_resolution_clock::now(); // --- stop
    if (N > 100000)
    {
      tree.printBT();
      tree.printInfo();
    }

    const auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
            .count();
    results.push_back(std::make_tuple(D, N, radius, time));

    std::cout << "DONE!";
    std::cout << "\nt: " << time << "ms";
    LINE_SEPARATOR;

    lsVTKWriter(mesh, lsFileFormatEnum::VTU, file_prefix + "TreeMesh.vtu")
        .apply();
  }
  LINE_SEPARATOR;
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
  test_2d.run(benchmark_sphere<2>);
  std::cout << "\n###################################################"
            << "\n######################## D = 3 ####################"
            << "\n###################################################"
            << std::endl;
  test_3d.run(benchmark_sphere<3>);
  LINE_SEPARATOR;

  std::cout << "------------- RESUMÉ -------------" << std::endl;

  check(test_2d, test_3d);

  if (all(test_2d.wasSuccess(), test_3d.wasSuccess()))
    std::cout << "Test " << argv[0] << ": SUCCESS!";
  else
    std::cout << "Test " << argv[0] << ": FAILURE!";

  return 0;
}
