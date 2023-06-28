#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <set>
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <cstdlib>
#include "lyra.hpp"
// ---------------------OpenMesh-----------------------------------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
//-----------------------------------------------------------------------------------------
using namespace OpenMesh;
using namespace Eigen;
struct MyTraits : public OpenMesh::DefaultTraits{
    VertexAttributes(OpenMesh::Attributes::Status);
};
typedef TriMesh_ArrayKernelT<MyTraits>  MyMesh;



class MeshWrap{
    public:
        MyMesh mesh;
    public:
        MeshWrap(const std::string& in, const std::string& out);
        ~MeshWrap();
        void initialize();
        void simplify(int n);
        void collapse_edge(MyMesh::EdgeHandle eh);
        void lock_boundary_edges();
        void recalculate_constraints_and_error();
        void write();
        void write_c(MyMesh::EdgeHandle eh);

        //for constraints calculation
        void get_constraints_and_error(MyMesh::EdgeHandle eh);
        bool is_alpha_compatible(MyMesh::EdgeHandle eh, Vector3d constraint);
        void add_constraint(MyMesh::EdgeHandle eh, Vector3d constraint, double b);
        void calc_remaining_constraints(MyMesh::EdgeHandle eh, MatrixXd Hessian, Vector3d c);
    private:
        int timestamp = 0;                  //Time measuring
        bool init = false;                  //Whether the error is initialized on all edges
        double alpha = 0.01745329251;       //angle, to which all planes are taken as coplanar (kind of)
        double lambda = 0.5;                //overall edge error equation constant variable
        std::string output_path;
        
        MyMesh::Normal              face_normal;
        MyMesh::Point               p, p0, p1;          //points for vertex coords storage
        MyMesh::Point               v_ideal_point;      //Handle point and Resulting vertex
        MyMesh::VertexIter          v_it, v_end;        //Vertex iterator
        MyMesh::EdgeIter            e_it, e_end;        //Edge iterator
        MyMesh::HalfedgeIter        he_it,he_end;       //Halfedge iterator
        MyMesh::VertexVertexIter    vv_it, vv_it_recalc;//Adjacent vertices of a vertex
        MyMesh::VertexEdgeIter      ve_it, ve_it_recalc;//Adjacent edge of a vertex
        MyMesh::VertexFaceIter      vf_it;              //Adjacent faces of a vertex
        MyMesh::FaceVertexIter      fv_it;              //Adjacent vertices of a face

        MyMesh::VertexHandle        vh1, vh2, vh;       //Vertex v1 to which v2 will collapse
        MyMesh::FaceHandle          fh1, fh2;           //Face handles
        MyMesh::EdgeHandle          eh, eh1, eh2;       //Edge handles
        MyMesh::HalfedgeHandle      heh;                //Halfedge handle
        std::vector<MyMesh::EdgeHandle>  eh_arr;        //Edge handle to store the error


        //Properties
        struct v_ideal_coords {Vector3d v = Vector3d::Zero(3);};
        struct b_vector {Vector3d b = Vector3d::Zero(3);};
        struct constraints {MatrixXd c = MatrixXd::Zero(3,3);};
        EPropHandleT<double>            e;          //Final edge error; will be sorted in ascending order
        EPropHandleT<v_ideal_coords>    v;          //Optimal coords for resulting vertex v 
        EPropHandleT<constraints>       a;          //Contraints for resulting vertex  v
        EPropHandleT<b_vector>          b;          //Contraint for resulting vertex  v, right side vector
        EPropHandleT<int>               n;          //number of constraints
        EPropHandleT<bool>              is_locked;  //applied to edges which have some adjacent vertices, that are locked, to conserve boundary. 

        Vector3d temp1 = Vector3d::Zero(3);         //temporary vector storage because cross function fails to resize
        Vector3d temp2 = Vector3d::Zero(3);         //temporary vector storage because cross function fails to resize
        Vector3d temp3 = Vector3d::Zero(3);         //temporary vector storage because cross function fails to resize


        

    public://time measurement
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        void time(std::chrono::time_point<std::chrono::high_resolution_clock> start);
};