#include "butterfly_scheme.h"

void butterflyScheme(Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
{
    auto newV = Eigen::MatrixXd{V.rows()*4,3}.setZero();
    auto newF = Eigen::MatrixXi{F.rows()*4,3}.setZero();
    
    //add old vertices
    for( int v = 0; v < V.rows(); ++v )
      newV.row(v) = V.row(v);

    //add new vertices
    int indexV = V.rows();
    int indexF = 0;
    std::map<std::pair<int,int>,int> m;
    for(int fi = 0; fi < F.rows(); ++fi) //iterate faces
    {
      auto f = F.row(fi);
      int point[3]={0,0,0};
      for(int i=0; i<3; i++) //iterate vertices in faces and add new vertex
      {
        if(!m.count(std::make_pair(f[i],f[(i+1)%3])))
        {
          point[i] = indexV++;
          newV.row(point[i]) = (V.row(f[i]) + V.row(f[(i+1)%3]))/2;
          m[std::make_pair(f[i],f[(i+1)%3])] = point[i];
          m[std::make_pair(f[(i+1)%3],f[i])] = point[i];
        }
        else
          point[i] = m[std::make_pair(f[i],f[(i+1)%3])];
      }
      //add four new faces
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(f[0], point[0], point[2]);
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(point[0], f[1], point[1]);
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(point[2], point[1], f[2]);
      newF.row(indexF++) = Eigen::Matrix<int,1,3>(point[2], point[0], point[1]);
    }
    newV.conservativeResize(indexV,3);
    
    trimesh::trimesh_t newMesh;
    buildMesh(newV, newF, newMesh);
    
    //update new vertices
    for(int vi = V.rows(); vi < newV.rows(); ++vi) 
    {
      std::vector< trimesh::index_t > neighs;
      std::vector< trimesh::index_t > old1,old2,old3,common(6);
      int v1,v2,v3,v4,v5,v6,v7,v8 = 0;
      //strange trick. The old vertices are always at index 0 and 3 of the neighbors
      newMesh.vertex_vertex_neighbors( vi, neighs );
      v1 = neighs.at(0);
      v2 = neighs.at(3);
      oldMesh.vertex_vertex_neighbors( v1, old1 );
      oldMesh.vertex_vertex_neighbors( v2, old2 );
      std::sort(old1.begin(), old1.end());
      std::sort(old2.begin(), old2.end());
      std::set_intersection(old1.begin(), old1.end(),old2.begin(), old2.end(),common.begin());
      v3 = common.at(0);
      v4 = common.at(1);
      oldMesh.vertex_vertex_neighbors( v3, old3 );
      std::sort(old1.begin(), old1.end());
      std::sort(old3.begin(), old3.end());
      std::set_intersection(old1.begin(), old1.end(),old3.begin(), old3.end(),common.begin());
      v5 = common.at(0) == v2 ? common.at(1) : common.at(0);
      std::sort(old2.begin(), old2.end());
      std::sort(old3.begin(), old3.end());
      std::set_intersection(old2.begin(), old2.end(),old3.begin(), old3.end(),common.begin());
      v6 = common.at(0) == v1 ? common.at(1) : common.at(0);
      oldMesh.vertex_vertex_neighbors( v4, old3 );
      std::sort(old1.begin(), old1.end());
      std::sort(old3.begin(), old3.end());
      std::set_intersection(old1.begin(), old1.end(),old3.begin(), old3.end(),common.begin());
      v7 = common.at(0) == v2 ? common.at(1) : common.at(0);
      std::sort(old2.begin(), old2.end());
      std::sort(old3.begin(), old3.end());
      std::set_intersection(old2.begin(), old2.end(),old3.begin(), old3.end(),common.begin());
      v8 = common.at(0) == v1 ? common.at(1) : common.at(0);

      //multliplying with 1/16 (instead of 1/8) provides MUCH better results, somehow
      Eigen::MatrixXd newVertexPos = (1.0/16.0)*(8.0*newV.row(v1)+8.0*newV.row(v2)+
        2.0*newV.row(v3)+2.0*newV.row(v4)-1.0*newV.row(v5)-1.0*newV.row(v6)-1.0*newV.row(v7)-1.0*newV.row(v8));
      newV.row(vi) = newVertexPos;
    }
    
    V = newV;
    F = newF;
    oldMesh = newMesh;
}

void butterflySchemeN(int N,Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
{
  for(int n=0;n<N;n++)
      butterflyScheme(V,F,oldMesh);
}