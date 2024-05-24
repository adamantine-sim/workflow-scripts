/* Copyright (c) 2023 - 2024, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/opencascade/manifold_lib.h>

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

// This is basically the same as dealii::OpenCASCADE closest_point
// but simplified for our needs. The interface also allows passing in
// a std::vector of faces instead of the whole shape.
dealii::Point<3> closest_point(const std::vector<TopoDS_Face> &faces,
                               const dealii::Point<3> &origin,
                               const double tolerance)
{
  gp_Pnt Pproj = dealii::OpenCASCADE::point(origin);

  double minDistance = std::numeric_limits<double>::max();
  gp_Pnt tmp_proj(0.0, 0.0, 0.0);

  for (const auto &face : faces)
  {
    // the projection function needs a surface, so we obtain the
    // surface upon which the face is defined
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

    ShapeAnalysis_Surface projector(SurfToProj);
    gp_Pnt2d proj_params =
        projector.ValueOfUV(dealii::OpenCASCADE::point(origin), tolerance);

    SurfToProj->D0(proj_params.X(), proj_params.Y(), tmp_proj);

    double distance = dealii::OpenCASCADE::point<3>(tmp_proj).distance(origin);
    if (distance < minDistance)
    {
      minDistance = distance;
      Pproj = tmp_proj;
    }
  }

  return dealii::OpenCASCADE::point<3>(Pproj);
}

// Try to move every boundary point to the surface described by the TopoDS_Faces
// passed in as argument. Optionally, ignore faces in z-direction.
void snap_to_iges(
    dealii::Triangulation<3> &tria, const std::vector<TopoDS_Face> &faces,
    dealii::OpenCASCADE::NormalProjectionManifold<3, 3> &projector,
    bool exclude_z_faces)
{

  // Store for every boundary point all its normal vectors. If all of these
  // point in the same direction, replace its x- and y-coordinates with the
  // respective coordinates of the closest point on the surface. Otherwise, the
  // point is a corner point and projecting it would collapse vertices. Instead
  // take the average of the old and new point.
  std::map<unsigned int, std::pair<std::reference_wrapper<dealii::Point<3>>,
                                   std::vector<dealii::Tensor<1, 3>>>>
      vertex_map;
  std::array<dealii::Tensor<1, 3>, dealii::GeometryInfo<3>::vertices_per_face>
      normal_at_vertex;
  for (const auto &cell : tria.active_cell_iterators())
  {
    for (const unsigned int i : cell->face_indices())
    {
      const auto &face = cell->face(i);
      if (face->at_boundary())
      {
        projector.get_normals_at_vertices(face, normal_at_vertex);
        for (unsigned j = 0; j < face->n_vertices(); ++j)
        {
          const unsigned int vertex_index = face->vertex_index(j);
          const auto &vertex_map_iterator = vertex_map.find(vertex_index);
          auto normal = normal_at_vertex[j] / normal_at_vertex[j].norm();

          if (!(exclude_z_faces &&
                std::abs(normal * dealii::Tensor<1, 3>{{0, 0, 1}}) > 0.1))
          {
            if (vertex_map_iterator == vertex_map.end())
            {
              std::pair<std::reference_wrapper<dealii::Point<3>>,
                        std::vector<dealii::Tensor<1, 3>>>
                  pair(face->vertex(j), {normal});
              vertex_map.emplace(vertex_index, pair);
            }
            else
            {
              std::get<1>(vertex_map_iterator->second).push_back(normal);
            }
          }
        }
      }
    }
  }

  for (const auto &boundary_vertex_iterator : vertex_map)
  {
    const auto &normals = std::get<1>(boundary_vertex_iterator.second);
    double minimum_product = 1;
    for (unsigned int i = 0; i < normals.size(); ++i)
    {
      for (unsigned int j = i + 1; j < normals.size(); ++j)
      {
        auto product = normals[i] * normals[j];
        minimum_product = std::min(product, minimum_product);
      }
    }
    auto &vertex = std::get<0>(boundary_vertex_iterator.second).get();
    auto proj = closest_point(faces, vertex, 1.e-10);
    // For tessellated meshes the minimal product between normal vectors can
    // only be -1,0, or 1.
    if (minimum_product > .5)
    {
      vertex(0) = proj(0);
      vertex(1) = proj(1);
    }
    else
    {
      vertex(0) = (vertex(0) + proj(0)) / 2;
      vertex(1) = (vertex(1) + proj(1)) / 2;
    }
  }
}

int main(int argc, char *argv[])
{
  if (argc < 4 || argc > 5)
  {
    std::cerr << "ERROR: The tool requires three runtime arguments for the "
                 "tessellated input vtk file, the IGES file, and the output "
                 "vtk file! Optionally, excluding z-facses can be requested. "
                 "However, the number of runtime arguments is "
              << argc - 1 << std::endl;
    std::abort();
  }
  std::cout << "Input VTK file: " << argv[1] << '\n'
            << "Input IGES file: " << argv[2] << '\n'
            << "Output VTK file: " << argv[3] << '\n';

  TopoDS_Shape surface = dealii::OpenCASCADE::read_IGES(argv[2]);

  double minima[3];
  double maxima[3];
  Bnd_Box B;
  BRepBndLib::Add(surface, B);
  B.Get(minima[0], minima[1], minima[2], maxima[0], maxima[1], maxima[2]);

  std::cout << "Bounding Box: (" << minima[0] << ',' << minima[1] << ','
            << minima[2] << "), (" << maxima[0] << ',' << maxima[1] << ','
            << maxima[2] << ")\n";

  dealii::OpenCASCADE::write_STL(surface, "debug.STL", 1.e-3);

  TopExp_Explorer exp;
  gp_Pnt tmp_proj;

  // Ignore TopoDS::Faces that are planes with z-normal. This might or might not
  // be a good choice depending on the input files. It should not be detrimental
  // if the face is indeed planar but might be problematic if we wrongly detect
  // it to be planar. Excluding those faces should mostly help with the edges at
  // the top and bottom of the mesh (assuming thet are planar).
  bool const exclude_z_faces = argc == 5 ? std::stoi(argv[4]) : false;

  std::vector<TopoDS_Face> faces;
  {
    for (exp.Init(surface, TopAbs_FACE); exp.More(); exp.Next())
    {
      TopoDS_Face face = TopoDS::Face(exp.Current());
      if (exclude_z_faces)
      {
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
        double aXmin, aYmin, aZmin, aXmax, aYmax, aZmax;

        Bnd_Box box;
        BRepBndLib::Add(face, box);
        box.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);

        // Obtain the four corners of the face and check if they are lying in
        // the same plane with z-normal
        double u1, u2, v1, v2;
        SurfToProj->Bounds(u1, u2, v1, v2);
        auto point_0 = dealii::OpenCASCADE::point<3>(SurfToProj->Value(u1, v1));
        auto point_1 = dealii::OpenCASCADE::point<3>(SurfToProj->Value(u1, v2));
        auto point_2 = dealii::OpenCASCADE::point<3>(SurfToProj->Value(u2, v1));
        auto point_3 = dealii::OpenCASCADE::point<3>(SurfToProj->Value(u2, v2));
        auto vector_0 = point_1 - point_0;
        auto vector_1 = point_2 - point_0;
        auto vector_2 = point_3 - point_0;
        auto normal = cross_product_3d(vector_0, vector_1);
        double deviation = normal / normal.norm() * vector_2 / vector_2.norm();
        double deviation_from_z =
            normal / normal.norm() * dealii::Tensor<1, 3>({0, 0, 1});

        // TODO find better tolerances to decide if a face is planar with
        // z-normal or not.
        if (std::abs(deviation) > .1 || std::abs(deviation_from_z) < .9)
          faces.push_back(face);
      }
      else
        faces.push_back(face);
    }
  }

  std::ifstream in;
  in.open(argv[1]);
  dealii::GridIn<3, 3> gi;
  dealii::Triangulation<3> tria;
  gi.attach_triangulation(tria);
  gi.read_vtk(in);

  dealii::BoundingBox vtk_bounding_box = tria.begin_active()->bounding_box();
  for (const auto &cell : tria.active_cell_iterators())
    vtk_bounding_box.merge_with(cell->bounding_box());
  std::cout << "Bounding Box (VTK): ("
            << vtk_bounding_box.get_boundary_points().first << ") ("
            << vtk_bounding_box.get_boundary_points().second << ")\n";

  // Compute shift and scaling to transform bouding boxes for debugging.
  for (int i = 0; i < 3; ++i)
  {
    double new_1 = vtk_bounding_box.get_boundary_points().first(i);
    double new_2 = vtk_bounding_box.get_boundary_points().second(i);
    double old_1 = minima[i];
    double old_2 = maxima[i];
    double shift = (old_2 * new_1 - new_2 * old_1) / (new_2 - new_1);
    double scale = new_1 / (old_1 + shift);
    std::cout << "shift[" << i << "]: " << shift << ", scale[" << i
              << "]: " << scale << '\n';
  }

  dealii::GridOut grid_out;
  {
    std::ofstream logfile("debug.vtk");
    grid_out.write_vtk(tria, logfile);
  }
  dealii::OpenCASCADE::NormalProjectionManifold<3, 3> normal_projector(surface);
  snap_to_iges(tria, faces, normal_projector, exclude_z_faces);

  const std::string filename = argv[3];
  std::ofstream logfile(filename);
  grid_out.write_vtk(tria, logfile);
}
