#include "loop_scheme.h"

void buildMesh(Eigen::MatrixXd& vertices, Eigen::MatrixXi& faces, trimesh::trimesh_t& mesh)
{
  std::vector< trimesh::triangle_t > triangles;
  std::vector< trimesh::edge_t > edges;
  int numV = vertices.rows();
  int numF = faces.rows();
  triangles.resize( numF );

  for (int i=0; i<numF; ++i)
  {
        triangles[i].v[0] = faces(i,0);
        triangles[i].v[1] = faces(i,1);
        triangles[i].v[2] = faces(i,2);
  }
  
  trimesh::unordered_edges_from_triangles( triangles.size(), &triangles[0], edges );
    
  mesh.build( numV, triangles.size(), &triangles[0], edges.size(), &edges[0] );
}

void loopScheme(Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
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
      std::vector< trimesh::index_t > old1,old2,common(6);
      int v1,v2,v3,v4 = 0;
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

      Eigen::MatrixXd newVertexPos = (1.0/8.0)*(3.0*newV.row(v1)+3.0*newV.row(v2)+newV.row(v3)+newV.row(v4));
      newV.row(vi) = newVertexPos;
    }

    //update old vertices
    for(int vi = 0; vi < V.rows(); ++vi) 
    {
      std::vector< trimesh::index_t > neighs;
      oldMesh.vertex_vertex_neighbors( vi, neighs );
      double n = neighs.size();
      double w = (64*n)/(40-pow(3+2*cos(2*PI/n),2))-n;
      auto sum = Eigen::MatrixXd{1,3}.setZero();
      for(int vk = 0; vk < neighs.size(); ++vk)
        sum += newV.row(neighs.at(vk));

      Eigen::MatrixXd newVertexPos = (1/(n+w))*(w * newV.row(vi) + sum);
      newV.row(vi) = newVertexPos;
    }
    V = newV;
    F = newF;
    oldMesh = newMesh;
}

void loopSchemeN(int N,Eigen::MatrixXd& V,Eigen::MatrixXi& F,trimesh::trimesh_t& oldMesh)
{
  for(int n=0;n<N;n++)
      loopScheme(V,F,oldMesh);

}