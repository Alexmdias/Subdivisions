#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include "trimesh.h"
#include "loop_scheme.h"
void sqrt3SchemeN(int N,Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh);