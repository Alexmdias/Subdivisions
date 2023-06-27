#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include "trimesh.h"
#define PI 3.141592653589793
void loopSchemeN(int N,Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh);
void buildMesh(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces, trimesh::trimesh_t& mesh);