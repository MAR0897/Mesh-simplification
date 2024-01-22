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
//#include <OpenMesh/Tools/Decimater/DecimaterT.hh>
//#include <OpenMesh/Tools/Decimater/ModQuadricT.hh>
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
        MeshWrap(const std::string& in, const std::string& out);    //Load mesh and activate mesh status, normals and properties
        ~MeshWrap();                                                //Output mesh into a designated file, deactivate mesh status, normals and properties
        void initialize();                                          //Initialize error and constraints on all edges of the mesh                    
        void simplify(int n);                                       //Decimate the mesh n-times (as specified in the input parameters)
        void collapse_edge(MyMesh::EdgeHandle eh);                  //Collapse one edge and delete it from the edge handle vector
        void recalculate_constraints_and_error();                   //Recalculates for all edges adjacent to vertices adjacent to the result collapsed vertex
        void lock_boundary_edges();                                 //Sets the "is_locked" edge property to true if edge has exactly 1 boundary vertex DOES NOT WORK
        void write();
        void write_c(MyMesh::EdgeHandle eh);

        void get_constraints_and_error(MyMesh::EdgeHandle eh);                                  //Computes the constraits (result vertex coords) and error of the edge
        bool is_alpha_compatible(MyMesh::EdgeHandle eh, Vector3d constraint);                   //Checks if to-be-added constraint is alpha compatible
        void add_constraint(MyMesh::EdgeHandle eh, Vector3d constraint, double b);              //Adds constraint to the system
        void calc_remaining_constraints(MyMesh::EdgeHandle eh, MatrixXd Hessian, Vector3d c);   //Calculates remaining constraints if there are less than 3

        double determinant3x3(const Matrix3d& mat);
        int prumer = 0;
    private:
        
        int timestamp = 0;                  //Time measuring
        bool init = false;                  //Whether the error is initialized on all edges
        double alpha = 0.01745329251;       //angle, to which all planes are taken as coplanar (kind of)
        double lambda = 0.5;                //overall edge error equation constant variable
        std::string output_path;
        
        MyMesh::Normal              face_normal;
        MyMesh::Point               p, p0, p1;          //points for vertex coords storage
        MyMesh::Point               v_ideal_point;      //Handle point for final collapsed vertex
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
        std::set<MyMesh::EdgeHandle> edge_handles_set; 

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

        //General
        Vector3d constraint = Vector3d::Zero(3);    //a - constraint
        double bside = 0;                           //b - right side number
        std::set<MyMesh::FaceHandle> face_handles;      //set of face handles (to avoid repeating calculations)
        std::set<MyMesh::VertexHandle> vertex_handles;  //set of vertex handles (to avoid repeating calculations)

        //Volume preservation
        double face_area;                           //face area storage
        double determinant;                         //face determinant storage
        Vector3d normal = Vector3d::Zero(3);        //face normal coords storage
        MatrixXd v_coords = MatrixXd::Zero(3,3);    //matrix for vertex coordinates to compute determinant

        //Volume, boundary optimization
        MatrixXd Hv = MatrixXd::Zero(3,3);          //Hessian for volume optimization
        MatrixXd Hb = MatrixXd::Zero(3,3);          //Hessian for boundary optimization
        MatrixXd Hs = MatrixXd::Zero(3,3);          //Hessian for triangle shape optimization
        Vector3d cv = Vector3d::Zero(3);            //for volume optimizaton
        Vector3d cb = Vector3d::Zero(3);            //for boundary optimization
        Vector3d cs = Vector3d::Zero(3);            //for triangle shape optimization
        double kv = 0;                              //constants in volume optimization
        double kb = 0;                              //constants in boundary optimization
        int i = 0;

        //Boundary preservation, optimization
        MatrixXd E1 = MatrixXd::Zero(3,3);          //e1 for every vertex
        MatrixXd E2 = MatrixXd::Zero(3,3);          //e2 for every vertex
        MatrixXd e1x = MatrixXd::Zero(3,3);         //e1x matrix for boundary optimization
        Vector3d e1 = Vector3d::Zero(3);            //summed E1
        Vector3d e2 = Vector3d::Zero(3);            //summed E2
        Vector3d e3 = Vector3d::Zero(3);            //cross product of e1 and e2
        Vector3d v0 = Vector3d::Zero(3);            //vertex coords storage
        Vector3d v1 = Vector3d::Zero(3);            //second vertex coords storage

        //Final error calculation
        double fv = 0;                              //volume objective functions (the resulting error)
        double fb = 0;                              //area objective functions (the resulting error)
        Vector3d V = Vector3d::Zero(3);             //final vertex

    public://time measurement
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        void time(std::chrono::time_point<std::chrono::high_resolution_clock> start);
        double cas = 0.0;
        double cas_collapse = 0.0;
        double cas_simp = 0.0;
        double cas_recalc = 0.0;
        double cas_inner = 0.0;
};