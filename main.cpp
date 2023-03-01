#include <iostream>
#include <fstream>                         
#include <vector>
#include <chrono>                          
#include <set>                              
#include <armadillo>                        //Matrices
#include <algorithm>                        //Sorting
#include <cmath>                            //Math library
// ---------------------OpenMesh-----------------------------------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
//-----------------------------------------------------------------------------------------
using namespace OpenMesh;
using namespace arma; //Petr Pavel reference
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
//-----------------------------------------------------------------------------------------
//Measures and prints time the code took
void time(std::chrono::time_point<std::chrono::high_resolution_clock> start);
    int timestamp = 0;
//Calculates Q error matrix for a vertex
void calc_Q(MyMesh &mesh, MyMesh::VertexHandle vh);
    double a,b,c,d,x,y,z;                       //Plane and point coords representation
    double make_abc_1_again;                    //Normalize a, b and c (a^2+b^2+c^2=1)
    double p[4];                                //Plane vector; p = [a,b,c,d]; plane = ax+by+cz+d=0
    MyMesh::Normal              normal;         //Normal of a face; allows for p vector calculation
    MyMesh::Point               point;          //Point
//Calculates ideal coords for vertex v
void calc_error_and_vcoords(MyMesh &mesh, MyMesh::EdgeHandle eh);
    MyMesh::HalfedgeHandle      heh;            //Function halfedge handle
    MyMesh::VertexHandle        vh1, vh2;       //Function vertex handle
    MyMesh::Point               midpoint;       //Midpoint in case optimal solution fails
//-----------------------------------------------------------------------------------------
    int i;
    bool init = false;                          //Was the error already calculated once? 
    double threshold = 0.002;                   //Threshold for connecting not connected vertices
    double mem[4];
    int rankA, rankAb;                          //Ranks of matrices; check if they are solvable
//-----------------------------------------------------------------------------------------
    MyMesh::Point               vvv;            //Handle point and Resulting vertex
    MyMesh::VertexIter          v_it, v_end;    //Vertex iterator
    MyMesh::EdgeIter            e_it, e_end;    //Edge iterator
    MyMesh::VertexVertexIter    vv_it;          //Adjacent vertices of a vertex
    MyMesh::VertexEdgeIter      ve_it;          //Adjacent edge of a vertex
    MyMesh::VertexFaceIter      vf_it;          //Adjacent faces of a vertex
    MyMesh::VertexHandle        v1, v;          //Vertex v1 to which v2 will collapse
    MyMesh::HalfedgeHandle      he;             //Edge handle
    MyMesh::EdgeHandle          eh;             //Edge handle
    std::vector<MyMesh::EdgeHandle>  eh_arr;    //Edge handle to store the error
    
    VPropHandleT<Error_matrix>      Q;          //Fundamental error matrix for each vertex
    EPropHandleT<double>            e;          //Final edge error; will be sorted in ascending order
    EPropHandleT<v_ideal_coords>    C;          //Optimal coords for vertex v

//Armadillo---------------------------------------------------------------------------------
    mat V(4, 4, fill::zeros);                   //Matrix for finding perfect spot for vertex v
    mat vektor(4,1, fill::zeros);               //Vector [0,0,0,1]
    mat v_coords(4,1, fill::zeros);             //Coordinates for the perfect spot for vertex v


//============================================================================================================================================================
//============================================================================================================================================================
int main(int argv, char **argc){
auto start = std::chrono::high_resolution_clock::now();

    vektor(3,0) = 1;
    std::ofstream myfile; //Temporary file check
    MyMesh mesh; if(!IO::read_mesh(mesh, "bunny.obj")) {std::cerr << "read error\n";exit(1);}       //Declare a triangular mesh object
    mesh.request_face_status();mesh.request_edge_status();mesh.request_vertex_status();             //Enable deleting objects
    mesh.request_face_normals();mesh.request_vertex_normals();                                      //Enable normals
    //Add properties
    mesh.add_property(Q);
    mesh.add_property(e);
    mesh.add_property(C);
    
//----------------------------------------------------------------------------------------------------------
    v_end=mesh.vertices_end();
    e_end=mesh.edges_end();

    //Compute initial Q error matrix for all vertices--------------------------------------------------------------------------
        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) calc_Q(mesh,*v_it);
    //Select all valid pairs-----------------------------------------------------------------------------------------------------
        //pro ted beru vsechny hrany, pozdeji bych mel i vrcholy ktere nejsou spojene ale jsou dost blizko k sobe
    //Compute the optimal vertex v-----------------------------------------------------------------------------------------------
        for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) calc_error_and_vcoords(mesh, *e_it);    
    //Sort eh_arr vector by smallest error----------------------------------------------------------------------------------------------
        std::sort(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh1, MyMesh::EdgeHandle eh2) {return mesh.property(e, eh1) < mesh.property(e, eh2);});
        
        //myfile.open ("output.txt");for (int m = 0; m<eh_arr.size(); m++) myfile << mesh.property(e, eh_arr[m]) << std::endl;myfile.close();
        init = true;
        time(start);

//INIT DONE        
//============================================================================================================================================================  
//Moving and deleting vertices
    for(int repeat = 0; repeat < atoi(argc[1]); repeat++){
    //Find mininal error
        eh = *std::min_element(eh_arr.begin(), eh_arr.end(),[&mesh](const MyMesh::EdgeHandle &eh_min1, const MyMesh::EdgeHandle &eh_min2) {
            return mesh.property(e, eh_min1) < mesh.property(e, eh_min2);});
    //Get its halfedge
        he = mesh.halfedge_handle(eh, 0);
    //Assign ideal coords    
        for (int j = 0; j<3; j++) vvv[j] = mesh.property(C, eh).v[j];
    //Collapse
        v1 = mesh.to_vertex_handle(he); mesh.set_point(v1, vvv); mesh.collapse(he);
    //Remove items marked as "deleted"   
        mesh.garbage_collection();
    //Erase from vector
        eh_arr.erase(std::remove_if(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh) {return !mesh.is_valid_handle(eh);}), eh_arr.end());

    //Recalculate error and ideal coords (prozatim na cele meshi)------------------------------------------------------------------------------------------
        if (mesh.is_valid_handle(v1)) {
            calc_Q(mesh, v1);       
            for (vv_it = mesh.vv_iter(v1); vv_it.is_valid(); ++vv_it) if(mesh.is_valid_handle(*vv_it)) {
                calc_Q(mesh, *vv_it);
                for (ve_it = mesh.ve_iter(*vv_it); ve_it.is_valid(); ++ve_it) if(!mesh.status(*ve_it).deleted()) {
                    calc_error_and_vcoords(mesh, *ve_it);
                }
            }
        }
   
    time(start);
    }          
//============================================================================================================================================================
//Output, deallocation

    if (!IO::write_mesh(mesh, "bunnyout.obj")) {std::cerr << "write error\n";exit(1);} 
    mesh.release_edge_status();mesh.release_vertex_status();mesh.release_face_status();mesh.release_face_normals();mesh.release_vertex_normals();
    mesh.remove_property(Q);
    mesh.remove_property(e);
    mesh.remove_property(C);

//Stop the timer
return 0;}
//============================================================================================================================================================
//============================================================================================================================================================
//FUNCTIONS 
//start je pouze jeden od nej se to odecita
void time(std::chrono::time_point<std::chrono::high_resolution_clock> start){
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "The code part #" << timestamp++ << " took " << duration.count() << " milliseconds to run." << std::endl;
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculate error matrix for some vertex
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
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculate ideal coords for vertex v
void calc_error_and_vcoords(MyMesh &mesh, MyMesh::EdgeHandle eh){
        for(int j = 0; j<4; j++) mem[j]=0.0;
    //get both vertices of the edge
        heh = mesh.halfedge_handle(eh, 0);
        vh1 = mesh.from_vertex_handle(heh);
        vh2 = mesh.to_vertex_handle(heh);
    //choose a perfect spot for vertex v (calculate ideal coords)
        for(int j = 0; j<3; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];
        V(3,0) = 0; V(3,1) = 0; V(3,2) = 0; V(3,3) = 1;
        rankA = rank(V);
        rankAb = rank(join_rows(V, vektor));
        if (rankA == rankAb) v_coords = solve(V,vektor/*,solve_opts::no_approx*/);
        else{
            midpoint = (mesh.point(vh1) + mesh.point(vh2)) / 2.0f;
            for (int j = 0; j<3; j++) v_coords(j,0) = midpoint[j];
        }
    //calculate error and write it to a vector of edge handles (using their properties)
        for(int j = 3; j<4; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];//plug the rest of the numbers to the error Q matrix (Q=Q1+Q2)
        if(!init) eh_arr.push_back(eh);
    //add error property
        mesh.property(e, eh) = 0.0;
        for(int j = 0; j<4; j++) for(int k = 0; k<4; k++) mem[j] += v_coords[k]*V(k,j);//matrix multiplication
        for(int j = 0; j<4; j++) mesh.property(e, eh) += mem[j]*v_coords[j];//matrix multiplication part 2
    //add ideal_v_coords property
        for(int j = 0; j<3; j++) mesh.property(C, eh).v[j] = v_coords(j,0);
}
