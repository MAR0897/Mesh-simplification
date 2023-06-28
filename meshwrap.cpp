#include "meshwrap.h"

        
MeshWrap::MeshWrap(const std::string& in, const std::string& out){
    //Read from file
    if (!IO::read_mesh(mesh, in)){std::cerr << "\n";exit(1);}
    output_path = out;
    //Request mesh status and normals, add properties
    mesh.request_face_status();mesh.request_edge_status();mesh.request_vertex_status();mesh.request_face_normals();mesh.request_vertex_normals();
    mesh.add_property(e); mesh.add_property(n); mesh.add_property(a); mesh.add_property(b); mesh.add_property(v);
    mesh.add_property(is_locked); 
    for (e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) mesh.property(is_locked, *e_it) = false;
}

MeshWrap::~MeshWrap(){
    //Write to file
    if (!IO::write_mesh(mesh, this->output_path)) {std::cerr << "write error\n";exit(1);} 
    //Release mesh status and normals, add properties
    mesh.release_edge_status();mesh.release_vertex_status();mesh.release_face_status();mesh.release_face_normals();mesh.release_vertex_normals();
    mesh.remove_property(e);
}

void MeshWrap::initialize(){
    e_end = mesh.edges_end();
    //for every edge get the constraints and the error
    for (e_it = mesh.edges_begin(); e_it != e_end; ++e_it) get_constraints_and_error(*e_it);
    init = true;
}

void MeshWrap::simplify(int n){
    for (int it = 0; it<n; it++){
        eh = *std::min_element(eh_arr.begin(), eh_arr.end(), [this](const MyMesh::EdgeHandle& eh_min1, const MyMesh::EdgeHandle& eh_min2) {
            return this->mesh.property(e, eh_min1) < this->mesh.property(e, eh_min2);});
        heh = mesh.halfedge_handle(eh, 0);
        vh1 = mesh.from_vertex_handle(heh);
        vh2 = mesh.to_vertex_handle(heh);
        if((mesh.is_boundary(vh1) xor mesh.is_boundary(vh2)) and (mesh.property(e, eh) == 0 and mesh.property(v, eh).v(0) == 0
        and mesh.property(v, eh).v(1) == 0 and mesh.property(v, eh).v(2) == 0)) {
            it--;
            eh_arr.erase(std::find(eh_arr.begin(), eh_arr.end(), eh));
        }
        else {
            collapse_edge(eh);
            recalculate_constraints_and_error();
        }
    }
}

void MeshWrap::collapse_edge(MyMesh::EdgeHandle eh){
    //Get its halfedge
    heh = mesh.halfedge_handle(eh, 0);
    //Assign ideal coords    
    for (int j = 0; j<3; j++) v_ideal_point[j] = mesh.property(v, eh).v[j];
    //Collapse
    vh1 = mesh.to_vertex_handle(heh);
    vh2 = mesh.from_vertex_handle(heh);
    mesh.set_point(vh1, v_ideal_point);mesh.collapse(heh);
    //Remove items marked as "deleted"   
    mesh.garbage_collection();
    //Erase from vector
    eh_arr.erase(std::remove_if(eh_arr.begin(), eh_arr.end(), [&](MyMesh::EdgeHandle eh) {return !mesh.is_valid_handle(eh);}), eh_arr.end());
}

void MeshWrap::recalculate_constraints_and_error() {
    /*std::set<MyMesh::EdgeHandle> edge_handles_set;
    if (mesh.is_valid_handle(vh1)) {  
        for (vv_it_recalc = mesh.vv_iter(vh1); vv_it_recalc.is_valid(); ++vv_it_recalc) if (mesh.is_valid_handle(*vv_it_recalc)) { 
            for (ve_it_recalc = mesh.ve_iter(*vv_it_recalc); ve_it_recalc.is_valid(); ++ve_it_recalc) if (!mesh.status(*ve_it_recalc).deleted()) { 
                edge_handles_set.insert(*ve_it_recalc);
            }
        }
    }
    for (const auto& edge_handle : edge_handles_set) get_constraints_and_error(edge_handle);
    edge_handles_set.clear();*/
    if (mesh.is_valid_handle(vh1)) {    
        for (vv_it_recalc = mesh.vv_iter(vh1); vv_it_recalc.is_valid(); ++vv_it_recalc) if(mesh.is_valid_handle(*vv_it_recalc)) { 
            for (ve_it_recalc = mesh.ve_iter(*vv_it_recalc); ve_it_recalc.is_valid(); ++ve_it_recalc) if(!mesh.status(*ve_it_recalc).deleted()) { 
                get_constraints_and_error(*ve_it_recalc);
            }
        }
    }
}


//==========================================================================================================================================================
void MeshWrap::time(std::chrono::time_point<std::chrono::high_resolution_clock> start){
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "The code part #" << timestamp++ << " took " << duration.count() << " milliseconds to run." << std::endl;
}

void MeshWrap::lock_boundary_edges(){
    for (he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end() ; ++he_it)
        if (mesh.is_boundary(*he_it) ) {
            vh1 = mesh.to_vertex_handle(*he_it);
            vh2 = mesh.from_vertex_handle(*he_it);
            mesh.status(vh1).set_locked(true);
            mesh.status(vh2).set_locked(true);
            for (ve_it = mesh.ve_iter(vh1); ve_it.is_valid(); ++ve_it) mesh.property(is_locked, *ve_it) = true;
            for (ve_it = mesh.ve_iter(vh2); ve_it.is_valid(); ++ve_it) mesh.property(is_locked, *ve_it) = true;
        }
    //for (e_it = mesh.edges_begin(); e_it != mesh.edges_end() ; ++e_it) if(mesh.is_boundary(*e_it)) if(mesh.is_valid_handle(*e_it)) mesh.status(*e_it).set_locked(true);
    
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
