#include "loop_scheme.h"
#include "butterfly_scheme.h"
#include "sqrt3_scheme.h"

int main(int argc, char *argv[])
{
  while(true)
  {
    int o,s,N = 1;
    std::cout << "Select Object (default : 1)\n";
    std::cout << 
      "cube: 1   armadillo: 2   bar: 3   bunny: 4   camel: 5   dino: 6    icosahedron: 7   sphere: 8    teddybear: 9\n";
    std::cin >> o;
    std::cout << "\nSelect N iteration (default : 1)\n";
    std::cin >> N;
    std::cout << "\nSelect Subdivision (default : 1)\n";
    std::cout << "loop: 1   butterfly: 2   sqrt3: 3\n";
    std::cin >> s;

    std::string str = "";
    switch(o)
    {
      case 1:
        str = "cube.obj";
        break;
      case 2:
        str = "armadillo_1k.off";
        break;
      case 3:
        str = "bar2.off";
        break;
      case 4:
        str = "bunny.obj";
        break;
      case 5:
        str = "camel_mc.obj";
        break;
      case 6:
        str = "dino.off";
        break;
      case 7:
        str = "icosahedron.obj";
        break;
      case 8:
        str = "sphere.obj";
        break;
      case 9:
        str = "teddybear.obj";
        break;
      default:
        str = "cube.obj";
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("../../../../input/"+str, V,F);
    trimesh::trimesh_t mesh;
    buildMesh(V,F,mesh);
    switch(s)
    {
      case 1:
        loopSchemeN(N,V,F,mesh);
        break;
      case 2:
        butterflySchemeN(N,V,F,mesh);
        break;
      case 3:
        sqrt3SchemeN(N,V,F,mesh);
        break;
      default:
        loopSchemeN(N,V,F,mesh);
    }

    // output the mesh
    igl::writeOBJ("../../../output/new"+str, V, F);
    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, mesh.get_faces());
    viewer.data().set_face_based(true);
    // add vertices highlights
    viewer.data().point_size = 20;
    viewer.data().add_points(V, Eigen::RowVector3d(1,0,0));
    // add vertices index
    for (int i=0; i<V.rows(); ++i)
        viewer.data().add_label(V.row(i)+Eigen::RowVector3d(0.005, 0.005, 0),std::to_string(i));
    viewer.data().show_custom_labels = true;
    // launch viewer
    viewer.launch();
  }
}