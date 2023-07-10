#include "meshwrap.h"

//Load mesh and activate mesh status, normals and properties        
MeshWrap::MeshWrap(const std::string& in, const std::string& out){
    if (!IO::read_mesh(mesh, in)){std::cerr << "\n";exit(1);}                               //Read the mesh from file
    output_path = out;                                                                      //Assign output_path variable
    mesh.request_face_status(); mesh.request_edge_status(); mesh.request_vertex_status();   //Request mesh status
    mesh.request_face_normals(); mesh.request_vertex_normals();                             //Request mesh normals
    mesh.add_property(e); mesh.add_property(n); mesh.add_property(a);                       //Add properties
    mesh.add_property(b); mesh.add_property(v); mesh.add_property(is_locked); 
    for (e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) mesh.property(is_locked, *e_it) = false;
}


//Output mesh into a designated file, deactivate mesh status, normals and properties
MeshWrap::~MeshWrap(){
    if (!IO::write_mesh(mesh, this->output_path)) {std::cerr << "write error\n";exit(1);}   //Write to file
    mesh.release_edge_status();mesh.release_vertex_status();mesh.release_face_status();     //Release mesh status 
    mesh.release_face_normals();mesh.release_vertex_normals();                              //Release mesh normals
    mesh.remove_property(e); mesh.remove_property(n); mesh.remove_property(a);              //Remove properties
    mesh.remove_property(b); mesh.remove_property(v); mesh.remove_property(is_locked);      
}


//Initialize error and constraints on all edges of the mesh 
void MeshWrap::initialize(){
    e_end = mesh.edges_end();                                                                   //Set the edge iterator end
    for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) get_constraints_and_error(*e_it);    //For every edge get the constraints and the error
    init = true;    //Set variable "init" to true, so that no more edges can get into the edge handle vector
}


//Decimate the mesh n-times (as specified in the input parameters)
void MeshWrap::simplify(int n){
    for (int it = 0; it<n; it++){                                                           //Iterate n-times (number_of_vertices parameter)
        eh = *std::min_element(eh_arr.begin(), eh_arr.end(), [this](const MyMesh::EdgeHandle& eh_min1, const MyMesh::EdgeHandle& eh_min2) {
            return this->mesh.property(e, eh_min1) < this->mesh.property(e, eh_min2);});    //Find the edge with minimal error
        heh = mesh.halfedge_handle(eh, 0);                                                  //Get its halfedge
        vh1 = mesh.from_vertex_handle(heh);
        vh2 = mesh.to_vertex_handle(heh);
        if((mesh.is_boundary(vh1) xor mesh.is_boundary(vh2)) and (mesh.property(e, eh) == 0 and mesh.property(v, eh).v(0) == 0
         and mesh.property(v, eh).v(1) == 0 and mesh.property(v, eh).v(2) == 0)) {          //If somehow the edge error and constraints are broken
            it--;                                                                           //  -decrement the number of vertices simplified
            eh_arr.erase(std::find(eh_arr.begin(), eh_arr.end(), eh));                      //  -erase the edge from edge handle vector
        }
        else {                                                                              //If all is OK
            collapse_edge(eh);                                                              //  -collapse the edge
            recalculate_constraints_and_error();                                            //  -recalculate the constrants and error where needed
        }
    }
}


//Collapse one edge and delete it from the edge handle vector
void MeshWrap::collapse_edge(MyMesh::EdgeHandle eh){
    heh = mesh.halfedge_handle(eh, 0);                                      //Get its halfedge
    for (int j = 0; j<3; j++) v_ideal_point[j] = mesh.property(v, eh).v[j]; //Assign ideal coords   
    vh1 = mesh.to_vertex_handle(heh);                                       
    vh2 = mesh.from_vertex_handle(heh);
    mesh.set_point(vh1, v_ideal_point);                                     //Move one of the vertices to the result vertex spot
    mesh.collapse(heh);                                                     //Collapse
    mesh.garbage_collection();                                              //Remove items marked as "deleted" 
    eh_arr.erase(std::remove_if(eh_arr.begin(), eh_arr.end(),               //Erase from vector of edge handles
    [&](MyMesh::EdgeHandle eh) {return !mesh.is_valid_handle(eh);}), eh_arr.end());
}


//Recalculates for all edges adjacent to vertices adjacent to the result collapsed vertex
void MeshWrap::recalculate_constraints_and_error() {
    std::set<MyMesh::EdgeHandle> edge_handles_set;                          //Create a set of edge handles and fill it with those, which need to recalculate their error
    if (mesh.is_valid_handle(vh1)) {  
        for (vv_it_recalc = mesh.vv_iter(vh1); vv_it_recalc.is_valid(); ++vv_it_recalc) if (mesh.is_valid_handle(*vv_it_recalc)) { 
            for (ve_it_recalc = mesh.ve_iter(*vv_it_recalc); ve_it_recalc.is_valid(); ++ve_it_recalc) if (!mesh.status(*ve_it_recalc).deleted()) { 
                edge_handles_set.insert(*ve_it_recalc);
            }
        }
    }
    for (const auto& edge_handle : edge_handles_set) get_constraints_and_error(edge_handle); //Recalculate the error
    edge_handles_set.clear();
}


//DOES NOT WORK Sets the "is_locked" edge property to true if edge has exactly 1 boundary vertex 
void MeshWrap::lock_boundary_edges(){
    for (he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); ++he_it) if(mesh.is_boundary(*he_it)) {
        vh1 = mesh.to_vertex_handle(*he_it);
        vh2 = mesh.from_vertex_handle(*he_it);
        mesh.status(vh1).set_locked(true);
        mesh.status(vh2).set_locked(true);
        for (ve_it = mesh.ve_iter(vh1); ve_it.is_valid(); ++ve_it) mesh.property(is_locked, *ve_it) = true;
        for (ve_it = mesh.ve_iter(vh2); ve_it.is_valid(); ++ve_it) mesh.property(is_locked, *ve_it) = true;
    }    
}

//==========================================================================================================================================================
void MeshWrap::time(std::chrono::time_point<std::chrono::high_resolution_clock> start){
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "The code part #" << timestamp++ << " took " << duration.count() << " milliseconds to run." << std::endl;
}
void MeshWrap::write(){
    std::sort(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh1, MyMesh::EdgeHandle eh2) {
    return mesh.property(e, eh1) < mesh.property(e, eh2);});
    std::ofstream outfile("output.txt");
    for (const auto& eh : eh_arr) outfile << mesh.property(e, eh) << "\t" << mesh.property(v, eh).v(0)
                     << "\t" << mesh.property(v, eh).v(1) << "\t" << mesh.property(v, eh).v(2) << "\n";
    outfile.close();
}
void MeshWrap::write_c(MyMesh::EdgeHandle eh){
    std::ofstream outfile2("output_c.txt", std::ios::app);
    outfile2 << mesh.property(b, eh).b << "\n" << mesh.property(a, eh).c << "\n\n\n";
    outfile2.close();
}

double MeshWrap::determinant3x3(const Matrix3d& mat){
    double det = 0.0;
    det += mat(0, 0) * (mat(1, 1) * mat(2, 2) - mat(2, 1) * mat(1, 2));
    det -= mat(0, 1) * (mat(1, 0) * mat(2, 2) - mat(2, 0) * mat(1, 2));
    det += mat(0, 2) * (mat(1, 0) * mat(2, 1) - mat(2, 0) * mat(1, 1));
    return det;
}
