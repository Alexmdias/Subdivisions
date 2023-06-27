#include "sqrt3_scheme.h"

void sqrt3Scheme(Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
{
    auto newV = Eigen::MatrixXd{V.rows()*4,3}.setZero();
    auto newF = Eigen::MatrixXi{F.rows()*3,3}.setZero();
    
    //add old vertices
    for( int v = 0; v < V.rows(); ++v )
      newV.row(v) = V.row(v);

    //add new vertices
    int indexV = V.rows();
    int indexF = 0;
    std::map<std::pair<int,int>,int> m;
    for(int fi = 0; fi < F.rows(); ++fi) //iterate faces
    {
      newV.row(indexV) = (1.0/3.0)*(V.row(F.row(fi)[0]) + V.row(F.row(fi)[1]) + V.row(F.row(fi)[2]));

      //add three new faces
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(indexV, F.row(fi)[0], F.row(fi)[1]);
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(indexV, F.row(fi)[1], F.row(fi)[2]);
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(indexV, F.row(fi)[2], F.row(fi)[0]);
      indexV++;
    }
    newV.conservativeResize(indexV,3);
    
    trimesh::trimesh_t newMesh;
    buildMesh(newV, newF, newMesh);
    
    //update old vertices
    for(int vi = 0; vi < V.rows(); ++vi) 
    {
      std::vector< trimesh::index_t > neighs;
      oldMesh.vertex_vertex_neighbors( vi, neighs );
      double n = neighs.size();
      auto sum = Eigen::MatrixXd{1,3}.setZero();
      for(int vk = 0; vk < neighs.size(); ++vk)
        sum += newV.row(neighs.at(vk));
      double a = (4-2*cos(2*PI/n))/9;
      Eigen::MatrixXd newVertexPos = (1.0-a)*newV.row(vi)+(a/n)*sum;
      newV.row(vi) = newVertexPos;
    }
    V = newV;
    F = newF;
    oldMesh = newMesh;
}

void sqrt3SchemeN(int N,Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
{
  for(int n=0;n<N;n++)
      sqrt3Scheme(V,F,oldMesh);
}