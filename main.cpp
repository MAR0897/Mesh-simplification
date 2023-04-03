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
    FaceAttributes(OpenMesh::Attributes::Color);
};
struct Error_matrix {
    double q[16];
};
struct v_ideal_coords {
    double v[3];
};
typedef TriMesh_ArrayKernelT<MyTraits>  MyMesh;
typedef Decimater::DecimaterT<MyMesh> MyDecimater;
typedef Decimater::ModQuadricT<MyMesh>::Handle MyHModQuadric;
//-----------------------------------------------------------------------------------------
//Measures and prints time the code took
void time(std::chrono::time_point<std::chrono::high_resolution_clock> start);
//-----------------------------------------------------------------------------------------
class MeshWrap{
    public:
        MyMesh mesh;
    public:
        MeshWrap(const std::string& in, const std::string& out){
            //Read from file
            if (!IO::read_mesh(mesh, in)){std::cerr << "\n";exit(1);}
            output_path = out;
            //Request mesh status and normals, add properties
            mesh.request_face_status();mesh.request_edge_status();mesh.request_vertex_status();mesh.request_face_normals();mesh.request_vertex_normals();
            mesh.add_property(Q);
            mesh.add_property(e);
            mesh.add_property(C);
            pomocny_vektor(3) = 1;
            v_end=mesh.vertices_end();
            e_end=mesh.edges_end();
        }
        ~MeshWrap(){
            //Write to file
            if (!IO::write_mesh(mesh, this->output_path)) {std::cerr << "write error\n";exit(1);} 
            //Release mesh status and normals, add properties
            mesh.release_edge_status();mesh.release_vertex_status();mesh.release_face_status();mesh.release_face_normals();mesh.release_vertex_normals();
            mesh.remove_property(Q);
            mesh.remove_property(e);
            mesh.remove_property(C);
        }
        //Calculates Q error matrix for a vertex
        void calc_Q(MyMesh &mesh, MyMesh::VertexHandle vh){
            //initialize all Q values to zero
            std::fill(std::begin(mesh.property(Q, vh).q), std::end(mesh.property(Q, vh).q), 0);
            for (vf_it = mesh.vf_iter(vh); vf_it.is_valid(); ++vf_it){
            //getting normal -> getting a, b, c
                normal=mesh.calc_face_normal(*vf_it);
                a = normal[0]; b = normal[1]; c = normal[2];
            //a^2+b^2+c^2=1
                make_abc_1_again = sqrt(1/(pow(a,2)+pow(b,2)+pow(c,2)));
                a *= make_abc_1_again; b *= make_abc_1_again; c *= make_abc_1_again;
            //getting coords of the vertex -> calculating d
                point=mesh.point(vh);
                x = point[0]; y = point[1]; z = point[2];
                d = -(a*x+b*y+c*z);
            //assinging to p vector -> computing Q matrix
                p[0] = a; p[1] = b; p[2] = c; p[3] = d;
                for(int j = 0; j<4; j++) for(int k = 0; k<4; k++) mesh.property(Q, vh).q[j*4+k] += p[j]*p[k];
            }
        }
        //Calculates ideal coords for vertex v
        void calc_error_and_vcoords(MyMesh &mesh, MyMesh::EdgeHandle eh){
                heh = mesh.halfedge_handle(eh, 0);
                vh1 = mesh.from_vertex_handle(heh);
                vh2 = mesh.to_vertex_handle(heh);
            //choose a perfect spot for vertex v (calculate ideal coords)
                for(int j = 0; j<3; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];
                V(3,0) = 0; V(3,1) = 0; V(3,2) = 0; V(3,3) = 1;
                if (1) v_coords = V.colPivHouseholderQr().solve(pomocny_vektor);
                else{
                    midpoint = (mesh.point(vh1) + mesh.point(vh2)) / 2.0f;
                    for (int j = 0; j<3; j++) v_coords(j) = midpoint[j];
                }
            //calculate error and write it to a vector of edge handles (using their properties)
                for(int j = 3; j<4; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];//plug the rest of the numbers to the error Q matrix (Q=Q1+Q2)
                if(!init) eh_arr.push_back(eh);
                mesh.property(e, eh) = v_coords.transpose()*(V*v_coords);//as_scalar(v_coords.t()*V*v_coords);
                if(mesh.is_boundary(vh1) || mesh.is_boundary(vh2)) mesh.property(e, eh)=DBL_MAX;
            //add ideal_v_coords property
                for(int j = 0; j<3; j++) mesh.property(C, eh).v[j] = v_coords(j);
        } 
        void get_error_and_ideal_coords(){
            //Compute initial Q error matrix for all vertices
            for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) calc_Q(mesh,*v_it);
            //Compute the optimal vertex v
            for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) calc_error_and_vcoords(mesh, *e_it); 
        }
        void collapse_one_edge(){
            //Find mininal error
                eh = *std::min_element(eh_arr.begin(), eh_arr.end(), [this](const MyMesh::EdgeHandle& eh_min1, const MyMesh::EdgeHandle& eh_min2) {
                    return this->mesh.property(e, eh_min1) < this->mesh.property(e, eh_min2);});
            //Get its halfedge
                he = mesh.halfedge_handle(eh, 0);
            //Assign ideal coords    
                for (int j = 0; j<3; j++) v_ideal_point[j] = mesh.property(C, eh).v[j];
            //Collapse
                v1 = mesh.to_vertex_handle(he);
                v2 = mesh.from_vertex_handle(he);
                mesh.set_point(v1, v_ideal_point);mesh.collapse(he);
            //Remove items marked as "deleted"   
                mesh.garbage_collection();
            //Erase from vector
                eh_arr.erase(std::remove_if(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh) {return !mesh.is_valid_handle(eh);}), eh_arr.end());
        }
        void recalculate_error_and_ideal_coords(){
            if (mesh.is_valid_handle(v1)) {
                calc_Q(mesh, v1);       
                for (vv_it = mesh.vv_iter(v1); vv_it.is_valid(); ++vv_it) if(mesh.is_valid_handle(*vv_it)) {
                    calc_Q(mesh, *vv_it);
                    for (ve_it = mesh.ve_iter(*vv_it); ve_it.is_valid(); ++ve_it) if(!mesh.status(*ve_it).deleted()) {
                        calc_error_and_vcoords(mesh, *ve_it);
                    }
                }
            }
        }
        void decimate(int vertices_to_decimate, int decimater_mode){
            switch(decimater_mode){
                case 1: {
                    MyDecimater decimater(mesh);
                    MyHModQuadric HQuadric;
                    decimater.add(HQuadric);
                    decimater.module(HQuadric).unset_max_err();
                    decimater.initialize();                    
                    decimater.decimate(vertices_to_decimate);
                    mesh.garbage_collection();
                    break;}
                case 2: {
                    get_error_and_ideal_coords(); 
                    init = true;
                    time(start);
                    for(int repeat = 0; repeat < vertices_to_decimate; repeat++){
                        collapse_one_edge();
                        recalculate_error_and_ideal_coords();
                        //time(start);
                    } 
                    break;}
            }
        }
        void lock_edges(){
            for (he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end() ; ++he_it)
                if (mesh.is_boundary(*he_it) ) {
                    mesh.status(mesh.to_vertex_handle(*he_it)).set_locked(true);
                    mesh.status(mesh.from_vertex_handle(*he_it)).set_locked(true);
                }
        }
    private:
        int timestamp = 0;                          //Time measuring
        bool init = false;                          //Will be set to TRUE when everything is initialized
        double a,b,c,d,x,y,z;                       //Plane and point coords representation
        double make_abc_1_again;                    //Normalize a, b and c (a^2+b^2+c^2=1)
        double p[4];                                //Plane vector; p = [a,b,c,d]; plane = ax+by+cz+d=0
        std::string output_path;
        MyMesh::Normal              normal;         //Normal of a face; allows for p vector calculation
        MyMesh::Point               point;          //Point
        MyMesh::HalfedgeHandle      heh;            //Function halfedge handle
        MyMesh::VertexHandle        vh1, vh2;       //Function vertex handle
        MyMesh::Point               midpoint;       //Midpoint in case optimal solution fails
        MyMesh::Point               v_ideal_point;  //Handle point and Resulting vertex
        MyMesh::VertexIter          v_it, v_end;    //Vertex iterator
        MyMesh::EdgeIter            e_it, e_end;    //Edge iterator
        MyMesh::HalfedgeIter        he_it,he_end;   //Halfedge iterator
        MyMesh::VertexVertexIter    vv_it;          //Adjacent vertices of a vertex
        MyMesh::VertexEdgeIter      ve_it;          //Adjacent edge of a vertex
        MyMesh::VertexFaceIter      vf_it;          //Adjacent faces of a vertex
        MyMesh::VertexHandle        v1, v2, v;          //Vertex v1 to which v2 will collapse
        MyMesh::HalfedgeHandle      he;             //Edge handle
        MyMesh::EdgeHandle          eh;             //Edge handle
        std::vector<MyMesh::EdgeHandle>  eh_arr;    //Edge handle to store the error

        VPropHandleT<Error_matrix>      Q;          //Fundamental error matrix for each vertex
        EPropHandleT<double>            e;          //Final edge error; will be sorted in ascending order
        EPropHandleT<v_ideal_coords>    C;          //Optimal coords for vertex v 

    //Eigen ---------------------------------------------------------------------------------
        MatrixXf V = MatrixXf::Zero(4,4);
        Vector4f pomocny_vektor = Vector4f::Zero(4);
        Vector4f v_coords = Vector4f::Zero(4);

        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    public:
        void time(std::chrono::time_point<std::chrono::high_resolution_clock> start){
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "The code part #" << timestamp++ << " took " << duration.count() << " milliseconds to run." << std::endl;
        }
};

//============================================================================================================================================================
//============================================================================================================================================================
int main(int argv, const char **argc){//./main -i input_file_path -o output_file_path -n vertices_to_decimate -d decimator_mode
 

    //Lyra arguments
    std::string input_path;
    std::string output_path;
    int number_of_vertices;
    int decimater_mode;
    bool show_help = false;
    //Parser
    auto cli
        = lyra::help(show_help)
        | lyra::opt(input_path, "input_path")
              ["-i"]["--input"]("Provide path for input .obj file")
        | lyra::opt(output_path, "output_path")
              ["-o"]["--output"]("Provide path for output .obj file")
        | lyra::opt(number_of_vertices, "number_of_vertices")
              ["-n"]["--number_of_vertices"]("How many vertices should the program decimate")
        | lyra::opt(decimater_mode, "decimater_mode")
              ["-d"]["--decimater"]("Choose the decimater you want to decimate your mesh with,\n1 = inbuilt OpenMesh quadric decimator\n2 = Garland decimation with locked edges");
    //Parse the program arguments, check, print
    auto result = cli.parse({argv, argc});
    if (show_help) {std::cout << cli << std::endl;return 0;}
    if (!result){std::cerr << "Error in command line: " << result.message() << std::endl;return 1;}
    std::cout << "Input path = " << input_path << ", Output path = " << output_path << ", Number of vertices = " << number_of_vertices << "\n";

    MeshWrap m(input_path, output_path);
    //m.lock_edges(); 
    m.decimate(number_of_vertices, decimater_mode);
          

return 0;}
