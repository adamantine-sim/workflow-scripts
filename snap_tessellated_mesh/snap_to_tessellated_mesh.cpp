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

#include <BRepBndLib.hxx>
#include <BRep_Tool.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>

#include <deal.II/opencascade/manifold_lib.h>

dealii::Point<3> my_closest_point(const std::vector<TopoDS_Face> &faces,
                                  const dealii::Point<3> &origin,
                                  const double tolerance) {
  gp_Pnt Pproj = dealii::OpenCASCADE::point(origin);

  double minDistance = std::numeric_limits<double>::max();
  gp_Pnt tmp_proj(0.0, 0.0, 0.0);

  for (const auto &face : faces) {
    // the projection function needs a surface, so we obtain the
    // surface upon which the face is defined
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);

    ShapeAnalysis_Surface projector(SurfToProj);
    gp_Pnt2d proj_params =
        projector.ValueOfUV(dealii::OpenCASCADE::point(origin), tolerance);

    SurfToProj->D0(proj_params.X(), proj_params.Y(), tmp_proj);

    double distance = dealii::OpenCASCADE::point<3>(tmp_proj).distance(origin);
    if (distance < minDistance) {
      minDistance = distance;
      Pproj = tmp_proj;
    }
  }

  return dealii::OpenCASCADE::point<3>(Pproj);
}

void snap_to_iges(
    dealii::Triangulation<3> &tria, const std::vector<TopoDS_Face> &faces,
    dealii::OpenCASCADE::NormalProjectionManifold<3, 3> &projector,
    bool exclude_z_faces) {

  std::map<unsigned int, std::pair<std::reference_wrapper<dealii::Point<3>>,
                                   std::vector<dealii::Tensor<1, 3>>>>
      vertex_map;
  std::array<dealii::Tensor<1, 3>, dealii::GeometryInfo<3>::vertices_per_face>
      normal_at_vertex;
  for (const auto &cell : tria.active_cell_iterators()) {
    for (const unsigned int i : cell->face_indices()) {
      const auto &face = cell->face(i);
      if (face->at_boundary()) {
        projector.get_normals_at_vertices(face, normal_at_vertex);
        for (unsigned j = 0; j < face->n_vertices(); ++j) {
          const unsigned int vertex_index = face->vertex_index(j);
          const auto &vertex_map_iterator = vertex_map.find(vertex_index);
          auto normal = normal_at_vertex[j] / normal_at_vertex[j].norm();

          if (!(exclude_z_faces &&
                std::abs(normal * dealii::Tensor<1, 3>{{0, 0, 1}}) < 0.1)) {
            if (vertex_map_iterator == vertex_map.end()) {
              std::pair<std::reference_wrapper<dealii::Point<3>>,
                        std::vector<dealii::Tensor<1, 3>>>
                  pair(face->vertex(j), {normal});
              vertex_map.emplace(vertex_index, pair);
            } else {
              std::get<1>(vertex_map_iterator->second).push_back(normal);
            }
          }
        }
      }
    }
  }

  for (const auto &boundary_vertex_iterator : vertex_map) {
    const auto &normals = std::get<1>(boundary_vertex_iterator.second);
    double minimum_product = 1;
    for (unsigned int i = 0; i < normals.size(); ++i) {
      for (unsigned int j = i + 1; j < normals.size(); ++j) {
        auto product = normals[i] * normals[j];
        minimum_product = std::min(product, minimum_product);
      }
    }
    auto &vertex = std::get<0>(boundary_vertex_iterator.second).get();
    auto proj = my_closest_point(faces, vertex, 1.e-10);
    if (minimum_product > .5) {
      vertex(0) = proj(0);
      vertex(1) = proj(1);
      // vertex(2) = proj(2); // hourglass
    } else {
      vertex(0) = (vertex(0) + proj(0)) / 2;
      vertex(1) = (vertex(1) + proj(1)) / 2;
      // vertex(2) = (vertex(2) + proj(2)) / 2; // hourglass
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr
        << "ERROR: The tool requires three runtime arguments for the "
           "tessellated input vtk file, the IGES file, and the output vtk file!"
        << std::endl;
    std::abort();
  }
  std::cout << "Input VTK file: " << argv[1] << '\n'
            << "Input IGES file: " << argv[2] << '\n'
            << "Output VTK file: " << argv[3] << '\n';

  TopoDS_Shape surface = dealii::OpenCASCADE::read_IGES(argv[2]);

  double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;
  Bnd_Box B;
  BRepBndLib::Add(surface, B);
  B.Get(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);

  std::cout << "Bounding Box: (" << Xmin << ',' << Ymin << ',' << Zmin << "), ("
            << Xmax << ',' << Ymax << ',' << Zmax << ")\n";

  TopExp_Explorer exp;
  gp_Pnt tmp_proj;

  bool exclude_z_faces = false;

  // We don't care about faces with z-normal so exclude them here
  std::vector<TopoDS_Face> faces;
  {
    for (exp.Init(surface, TopAbs_FACE); exp.More(); exp.Next()) {
      TopoDS_Face face = TopoDS::Face(exp.Current());
      if (exclude_z_faces) {
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
        double aXmin, aYmin, aZmin, aXmax, aYmax, aZmax;

        Bnd_Box box;
        BRepBndLib::Add(face, box);
        box.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);

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

        if (std::abs(deviation) > .1 || std::abs(deviation_from_z) < .9)
          faces.push_back(face);
      } else
        faces.push_back(face);
    }
  }

  std::ifstream in;
  in.open(argv[1]);
  dealii::GridIn<3, 3> gi;
  dealii::Triangulation<3> tria;
  gi.attach_triangulation(tria);
  gi.read_vtk(in);
  // For Hourglass
  /*dealii::GridTools::transform ([](const dealii::Point<3> &p) ->
     dealii::Point<3>
                    {
                      dealii::Point<3> q = p;
                      std::swap(q[1], q[2]);
                      q[0] -= 0.07;
                      q[2] -= 0.07;
                      return q;
                    },
                    tria);*/
  // For curved wall
  dealii::GridTools::scale(.01, tria);

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
