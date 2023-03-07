#include <iostream>
#include <fstream>                          //Files
#include <vector>
#include <chrono>                           //Time
#include <set>                              //List
#include <armadillo>                        //Matrices
#include <algorithm>                        //Sorting
#include <cmath>                            //Normal math library
// ---------------------OpenMesh------------------------------------------------------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <OpenMesh/Core/Mesh/Attributes.hh>
//#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
//------------------------------------------------------------------------------------------------------------
using namespace OpenMesh;
using namespace arma;
struct MyTraits : public OpenMesh::DefaultTraits{
};
typedef TriMesh_ArrayKernelT<MyTraits>  MyMesh;
struct Error_matrix {
    double q[16];
};
struct v_ideal_coords {
    double v[3];
};
typedef struct Unconnected_vertices{
    MyMesh::VertexHandle uv1, uv2;
    double error;
    double v[3];
}Unconnected_vertices;
//------------------------------------------------------------------------------------------------------------
//Measures and prints time the code took
void time(std::chrono::time_point<std::chrono::high_resolution_clock> start);
    int timestamp = 0;
//Calculates Q error matrix for a vertex
void calc_Q(MyMesh &mesh, MyMesh::VertexHandle vh);
    double a,b,c,d,x,y,z;                       //Plane and point coords representation
    double make_abc_1_again;                    //Normalize a, b and c (a^2+b^2+c^2=1)
    double p[4];                                //Plane vector; p = [a,b,c,d]; plane = ax+by+cz+d=0
//Calculates ideal coords for vertex v
void calc_error_and_vcoords(MyMesh &mesh, MyMesh::EdgeHandle eh);    
//Distance between two vertices
double distance(MyMesh &mesh, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2);
void collapse_edge_only(MyMesh &mesh, MyMesh::HalfedgeHandle hh);
void uv_calc_error_and_vcoords(MyMesh &mesh, Unconnected_vertices &u_vertex, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2);
//------------------------------------------------------------------------------------------------------------
    bool init;                                  //Will be set to TRUE after INIT is done
    double threshold;                           //Threshold for connecting not connected vertices
    double mem[4];                              //Matrix multiplicaiton memory
    int rankA, rankAb;                          //Ranks of matrices; check if they are solvable
//OpenMesh----------------------------------------------------------------------------------------------------
    MyMesh::Point               v_ideal_point, point, midpoint;   //Handle point and Resulting vertex
    MyMesh::Normal              normal;                 //Normal for calculating a, b, c
    MyMesh::VertexIter          v_it, v_it2, v_end;     //Vertex iterator
    MyMesh::EdgeIter            e_it, e_end;            //Edge iterator
    MyMesh::VertexVertexIter    vv_it;                  //Adjacent vertices of a vertex
    MyMesh::VertexEdgeIter      ve_it;                  //Adjacent edges of a vertex
    MyMesh::VertexFaceIter      vf_it;                  //Adjacent faces of a vertex
    MyMesh::VertexHandle        vh1, vh2, vh;           //Vertices v1 and v2 which will collapse to only one vertex v
    MyMesh::HalfedgeHandle      heh;                    //Halfedges for v1 and v2
    MyMesh::EdgeHandle          eh;                     //Edge handle
    std::vector<MyMesh::EdgeHandle>  eh_arr;            //Edge handle to store the error
    VPropHandleT<Error_matrix>      Q;                  //Fundamental error matrix for each vertex
    EPropHandleT<double>            e;                  //Final edge error; will be sorted in ascending order
    EPropHandleT<v_ideal_coords>    C;                  //Optimal coords for vertex v
//Armadillo----------------------------------------------------------------------------------------------------
    mat V(4, 4, fill::zeros);                           //Matrix for finding perfect spot for vertex v
    mat pomocny_vektor(4, 1, fill::zeros);              //Vector [0,0,0,1]
    mat v_coords(4,1, fill::zeros);                     //Coordinates for the perfect spot for vertex v


//============================================================================================================================================================
//============================================================================================================================================================
int main(int argv, char **argc){//  ./main [number of vertices to delete] [threshold]
auto start = std::chrono::high_resolution_clock::now();

    threshold = atof(argc[2]);
    pomocny_vektor(3,0) = 1;
    std::ofstream myfile; //Temporary file check
    std::vector<Unconnected_vertices> uv;
    MyMesh mesh; if(!IO::read_mesh(mesh, "bunny.obj")) {std::cerr << "read error\n";exit(1);}       //Declare a triangular mesh object
    mesh.request_face_status();mesh.request_edge_status();mesh.request_vertex_status();             //Enable deleting objects
    mesh.request_face_normals();mesh.request_vertex_normals();                                      //Enable normals
    //Add properties
    mesh.add_property(Q);
    mesh.add_property(e);
    mesh.add_property(C);
    
    int i = 0;
//------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Compute initial Q error matrix for all vertices
        v_end=mesh.vertices_end();
        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) calc_Q(mesh,*v_it);
    //Select all valid pairs (those not connected)
        for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it){
            for (v_it2 = v_it+1; v_it2 != v_end; ++v_it2){
                if (mesh.find_halfedge(*v_it, *v_it2) == MyMesh::HalfedgeHandle()) {
                    if (threshold >= distance(mesh, *v_it, *v_it2) != 0) {
                        Unconnected_vertices u_vertex;
                        uv_calc_error_and_vcoords(mesh, u_vertex, *v_it, *v_it2);
                        //heh = mesh.new_edge(*v_it, *v_it2);
                        //U kralika threshold 0.01 vytvori celkem ~18000 edgu (z toho ~7500 je real)
                        uv.push_back(u_vertex);
                    }
                }
            }
        }  
    //Compute the optimal vertex v
        e_end=mesh.edges_end();
        for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) calc_error_and_vcoords(mesh, *e_it);
    //Sort eh_arr vector by smallest error
        std::sort(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh1, MyMesh::EdgeHandle eh2) {return mesh.property(e, eh1) < mesh.property(e, eh2);});
        init = true;


time(start);
//INIT DONE        
//============================================================================================================================================================  
//Merging two vertices into one argc[1] times
    for(int repeat = 0; repeat < atoi(argc[1]); repeat++){
    //Moving vertex and collapsing edge
        //Choose the edge with smallest error
        auto uv_it = std::min_element(uv.begin(), uv.end(), [](const Unconnected_vertices &a, const Unconnected_vertices &b) {return a.error < b.error;});
        eh = *std::min_element(eh_arr.begin(), eh_arr.end(), [&mesh](const MyMesh::EdgeHandle &eh_min1, const MyMesh::EdgeHandle &eh_min2) {
                    return mesh.property(e, eh_min1) < mesh.property(e, eh_min2);});

//------------------------------------------------------------------------------------------------------------------------------------------------------------
    //CONNECTED VERTICES
        if(mesh.property(e, eh) < uv_it->error){//Compare connected and unconnected vertices
            heh = mesh.halfedge_handle(eh, 0);
            //Assign ideal coords     
            for (int j = 0; j<3; j++) v_ideal_point[j] = mesh.property(C, eh).v[j];
            //Move one vertex   
            vh = mesh.to_vertex_handle(heh); mesh.set_point(vh, v_ideal_point);
            mesh.collapse(heh);
            //Remove items marked as "deleted"
            mesh.garbage_collection();
            //Erase deleted edges from eh_arr vector
            eh_arr.erase(std::remove_if(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh) {return !mesh.is_valid_handle(eh);}), eh_arr.end());
        //Recalculate error and ideal coords
            if (mesh.is_valid_handle(vh)) {
                calc_Q(mesh, vh);      
                for (vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) if(mesh.is_valid_handle(*vv_it)) {
                    calc_Q(mesh, *vv_it);
                    for (ve_it = mesh.ve_iter(*vv_it); ve_it.is_valid(); ++ve_it) if(!mesh.status(*ve_it).deleted()) {
                        calc_error_and_vcoords(mesh, *ve_it);}}}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
    //UNCONNECTED VERTICES
        }else{
            heh = mesh.new_edge(uv_it->uv1, uv_it->uv2);          
/*MOZNA ZBYTECNE*/eh = mesh.edge_handle(heh);
/*MOZNA ZBYTECNE*/heh = mesh.halfedge_handle(eh, 0);
            //Assign ideal coords 
            for (int j = 0; j<3; j++) v_ideal_point[j] = uv_it->v[j];           
            //Erase uv element containing second vertex
            uv.erase(std::remove_if(uv.begin(), uv.end(),[uv_it](const Unconnected_vertices &uv_elem) {return uv_elem.uv2 == uv_it->uv2;}),uv.end());
            uv.erase(std::remove_if(uv.begin(), uv.end(),[uv_it](const Unconnected_vertices &uv_elem) {return uv_elem.uv1 == uv_it->uv1;}),uv.end());
            //Move one vertex   
            vh = mesh.to_vertex_handle(heh); mesh.set_point(vh, v_ideal_point);
            if(mesh.next_halfedge_handle(heh).is_valid()) collapse_edge_only(mesh, heh);
            mesh.delete_edge(eh);
            //Remove items marked as "deleted"
            mesh.garbage_collection();
        //Recalculate error and ideal coords
            if (mesh.is_valid_handle(vh)) {
                calc_Q(mesh, vh);     
                for (vv_it = mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it) if(mesh.is_valid_handle(*vv_it)) {
                    calc_Q(mesh, *vv_it);
                    for (ve_it = mesh.ve_iter(*vv_it); ve_it.is_valid(); ++ve_it) if(!mesh.status(*ve_it).deleted()) {
                        calc_error_and_vcoords(mesh, *ve_it);
                    }
                }
                v_end = mesh.vertices_end();
                for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it){
                    if (mesh.find_halfedge(*v_it, vh) == MyMesh::HalfedgeHandle()) {
                        if (threshold >= distance(mesh, *v_it, vh) != 0) {
                            Unconnected_vertices u_vertex;
                            uv_calc_error_and_vcoords(mesh, u_vertex, *v_it, vh);
                            uv.push_back(u_vertex);
                        }
                    }
                }
            }
        }
    //std::cout<<"ahojda"<<std::endl;            
    time(start);
    }         
//===========================================================================================================================================================
//Output, deallocation 
std::cout<<"ahojda"<<std::endl;
    if (!IO::write_mesh(mesh, "bunnyout.obj")) {std::cerr << "write error\n";exit(1);} 
    mesh.release_edge_status();mesh.release_vertex_status();mesh.release_face_status();mesh.release_face_normals();mesh.release_vertex_normals();
    mesh.remove_property(Q);
    mesh.remove_property(e);
    mesh.remove_property(C);

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
//distance between two vertices
double distance(MyMesh &mesh, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2){
    MyMesh::Point p1 = mesh.point(v1);
    MyMesh::Point p2 = mesh.point(v2);
    return sqrt(pow(p1[0]-p2[0], 2) + pow(p1[1]-p2[1], 2) + pow(p1[2]-p2[2], 2));
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculate error matrix for some vertex
void calc_Q(MyMesh &mesh, MyMesh::VertexHandle vh){
    //Initialize all Q values to zero
    std::fill(std::begin(mesh.property(Q, vh).q), std::end(mesh.property(Q, vh).q), 0);
    for (vf_it = mesh.vf_iter(vh); vf_it.is_valid(); ++vf_it){
    //Getting normal -> getting a, b, c and then a^2+b^2+c^2=1
        normal=mesh.calc_face_normal(*vf_it);
        a = normal[0]; b = normal[1]; c = normal[2];
        make_abc_1_again = sqrt(1/(pow(a,2)+pow(b,2)+pow(c,2)));
        a *= make_abc_1_again; b *= make_abc_1_again; c *= make_abc_1_again;
    //Getting coords of the vertex -> calculating d
        point=mesh.point(vh);
        x = point[0]; y = point[1]; z = point[2];
        d = -(a*x+b*y+c*z);
    //Assinging to p vector -> computing Q matrix
        p[0] = a; p[1] = b; p[2] = c; p[3] = d;
        for(int j = 0; j<4; j++) for(int k = 0; k<4; k++) mesh.property(Q, vh).q[j*4+k] += p[j]*p[k];
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculate ideal coords for vertex v
void calc_error_and_vcoords(MyMesh &mesh, MyMesh::EdgeHandle eh){
        for(int j = 0; j<4; j++) mem[j]=0.0;
        if(!init) eh_arr.push_back(eh);
    //Get both vertices of the edge
        heh = mesh.halfedge_handle(eh, 0);
        vh1 = mesh.from_vertex_handle(heh);
        vh2 = mesh.to_vertex_handle(heh);
    //Calculate ideal coords and add it as an edge property
        for(int j = 0; j<3; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];
        V(3,0) = 0; V(3,1) = 0; V(3,2) = 0; V(3,3) = 1;
        rankA = rank(V);
        rankAb = rank(join_rows(V, pomocny_vektor));
        if (rankA == rankAb) v_coords = solve(V,pomocny_vektor/*,solve_opts::no_approx*/);
        else{
            midpoint = (mesh.point(vh1) + mesh.point(vh2)) / 2.0f;
            for (int j = 0; j<3; j++) v_coords(j,0) = midpoint[j];}
        for(int j = 0; j<3; j++) mesh.property(C, eh).v[j] = v_coords(j,0);
    //Calculate error and add it as an edge property
        for(int j = 3; j<4; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, vh1).q[j*4+k] + mesh.property(Q, vh2).q[j*4+k];//plug the rest of the numbers to the error Q matrix (Q=Q1+Q2)
        mesh.property(e, eh) = 0.0;
        for(int j = 0; j<4; j++) for(int k = 0; k<4; k++) mem[j] += v_coords[k]*V(k,j);//matrix multiplication
        for(int j = 0; j<4; j++) mesh.property(e, eh) += mem[j]*v_coords[j];//matrix multiplication part 2    
}
//The same function but for unconnected vertices
void uv_calc_error_and_vcoords(MyMesh &mesh, Unconnected_vertices &u_vertex, MyMesh::VertexHandle v1, MyMesh::VertexHandle v2){
        for(int j = 0; j<4; j++) mem[j]=0.0;
        int temp = 0;
        u_vertex.uv1 = v1;
        u_vertex.uv2 = v2;
    //Calculate ideal coords and add it as an edge property
        for(int j = 0; j<3; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, v1).q[j*4+k] + mesh.property(Q, v2).q[j*4+k];
        V(3,0) = 0; V(3,1) = 0; V(3,2) = 0; V(3,3) = 1;
        rankA = rank(V);
        rankAb = rank(join_rows(V, pomocny_vektor));
        if (rankA == rankAb) v_coords = solve(V,pomocny_vektor/*,solve_opts::no_approx*/);
        else{
            midpoint = (mesh.point(v1) + mesh.point(v2)) / 2.0f;
            for (int j = 0; j<3; j++) v_coords(j,0) = midpoint[j];}
        for(int j = 0; j<3; j++) u_vertex.v[j] = v_coords(j,0);
    //Calculate error and add it as an edge property
        for(int j = 3; j<4; j++) for(int k = 0; k<4; k++) V(j,k) = mesh.property(Q, v1).q[j*4+k] + mesh.property(Q, v2).q[j*4+k];//plug the rest of the numbers to the error Q matrix (Q=Q1+Q2)
        for(int j = 0; j<4; j++) for(int k = 0; k<4; k++) mem[j] += v_coords[k]*V(k,j);//matrix multiplication
        for(int j = 0; j<4; j++) temp += mem[j]*v_coords[j];//matrix multiplication part 2 
        u_vertex.error = temp;
        
}
void collapse_edge_only(MyMesh &mesh, MyMesh::HalfedgeHandle hh){
    MyMesh::HalfedgeHandle  h  = hh;
    MyMesh::HalfedgeHandle  hn = mesh.next_halfedge_handle(h);
    MyMesh::HalfedgeHandle  hp = mesh.prev_halfedge_handle(h);
    MyMesh::HalfedgeHandle  o  = mesh.opposite_halfedge_handle(h);
    MyMesh::HalfedgeHandle  on = mesh.next_halfedge_handle(o);
    MyMesh::HalfedgeHandle  op = mesh.prev_halfedge_handle(o);
    MyMesh::FaceHandle      fh = mesh.face_handle(h);
    MyMesh::FaceHandle      fo = mesh.face_handle(o);
    MyMesh::VertexHandle    vh = mesh.to_vertex_handle(h);
    MyMesh::VertexHandle    vo = mesh.to_vertex_handle(o);
    // halfedge -> vertex
    for (MyMesh::VertexIHalfedgeIter vih_it(mesh.vih_iter(vo)); vih_it.is_valid(); ++vih_it)
        mesh.set_vertex_handle(*vih_it, vh);
    // halfedge -> halfedge
    mesh.set_next_halfedge_handle(hp, hn);
    mesh.set_next_halfedge_handle(op, on);
    // face -> halfedge
    if (fh.is_valid())  mesh.set_halfedge_handle(fh, hn);
    if (fo.is_valid())  mesh.set_halfedge_handle(fo, on);
    // vertex -> halfedge
    if (mesh.halfedge_handle(vh) == o)  mesh.set_halfedge_handle(vh, hn);
    mesh.adjust_outgoing_halfedge(vh);
    mesh.set_isolated(vo);
    // delete stuff
    mesh.status(mesh.edge_handle(h)).set_deleted(false);
    mesh.status(vo).set_deleted(false);
    if (mesh.has_halfedge_status()){
        mesh.status(h).set_deleted(false);
        mesh.status(o).set_deleted(false);}
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------
/*myfile.open ("vcoords.txt");
        for (int m = 0; m<eh_arr.size(); m++) myfile << mesh.property(C, eh_arr[m]).v[0] <<"\t" <<mesh.property(C, eh_arr[m]).v[1]<<"\t"<< 
        mesh.property(C, eh_arr[m]).v[2]<< std::endl;
        myfile.close();*/

/*myfile.open ("output.txt");for (int m = 0; m<eh_arr.size(); m++) myfile << mesh.property(e, eh_arr[m]) << std::endl;myfile.close();*/

    /*myfile.open ("output2.txt");for (int m = 0; m<eh_arr.size(); m++) {
                if(!(mesh.face_handle(mesh.halfedge_handle(eh_arr[m], 0)).is_valid() || mesh.face_handle(mesh.halfedge_handle(eh_arr[m], 1)).is_valid())) {
                    myfile << mesh.property(e, eh_arr[m]) << "\tisolated";}
                else if(mesh.is_boundary(eh_arr[m])) myfile << mesh.property(e, eh_arr[m]) << "\tboundary";
                else myfile << mesh.property(e, eh_arr[m]) << "\tnormaaal";    
                if (mesh.status(eh_arr[m]).deleted()) {myfile << "\t" << m << "\tdeleted" << std::endl;} 
                else {myfile << "\t" << m << "\tvalid" << std::endl;}} myfile.close();*/