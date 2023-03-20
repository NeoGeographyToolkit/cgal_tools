#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>

#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>

#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <cstdlib>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Sphere_3 Sphere_3;
typedef CGAL::Point_set_3<Point_3, Vector_3> Point_set;

// TODO(oalexan1): This needs more testing.
// Usage: poisson_remesh input.ply output.ply

int main(int argc, char*argv[]) {

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << "<mode> input.ply output.ply"
              << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_file = argv[1];
  const char* output_file = argv[2]; 

  std::ifstream stream (input_file, std::ios_base::binary);
  if (!stream) {
    std::cerr << "Error: cannot read file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  Point_set points;
  stream >> points;

  std::cout << "Read " << points.size () << " point(s)" << std::endl;
  if (points.empty())
    return EXIT_FAILURE;

  // Outlier removal
  CGAL::remove_outliers<CGAL::Sequential_tag>
    (points,
     24, // Number of neighbors considered for evaluation
     points.parameters().threshold_percent (5.0)); // Percentage of points to remove
  
  std::cout << points.number_of_removed_points()
            << " point(s) are outliers." << std::endl;
  
  // Applying point set processing algorithm to a CGAL::Point_set_3
  // object does not erase the points from memory but place them in
  // the garbage of the object: memory can be freeed by the user.
  points.collect_garbage();

  // Simplification

  // Compute average spacing using neighborhood of 6 points
  double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag> (points, 6);

  // Simplify using a grid of size 2 * average spacing
  CGAL::grid_simplify_point_set(points, 2. * spacing);

  std::cout << points.number_of_removed_points()
            << " point(s) removed after simplification." << std::endl;

  points.collect_garbage();

  // Smoothing
  CGAL::jet_smooth_point_set<CGAL::Sequential_tag> (points, 24);

  CGAL::Surface_mesh<Point_3> output_mesh;

  std::cout << "Estimate normals" << std::endl;
  CGAL::jet_estimate_normals<CGAL::Sequential_tag>
    (points, 24); // Use 24 neighbors
  
  // Orientation of normals, returns iterator to first unoriented point
  typename Point_set::iterator unoriented_points_begin =
    CGAL::mst_orient_normals(points, 24); // Use 24 neighbors
  
  std::cout << "Clean up normals" << std::endl;
  points.remove (unoriented_points_begin, points.end());
  
  // Poisson reconstruction
  std::cout << "Poisson surface reconstruction" << std::endl;  
  CGAL::poisson_surface_reconstruction_delaunay
    (points.begin(), points.end(),
     points.point_map(), points.normal_map(),
     output_mesh, spacing);
  
  std::cout << "Writing output mesh: " << output_file << std::endl;
  CGAL::IO::write_PLY(output_file, output_mesh,
                      CGAL::parameters::stream_precision(17).use_binary_mode(false));
  
  
  return EXIT_SUCCESS;
}
