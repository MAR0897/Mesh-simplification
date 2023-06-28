#include "meshwrap.h"


void MeshWrap::get_constraints_and_error(MyMesh::EdgeHandle eh){

    if(!mesh.property(is_locked, eh)){

        mesh.property(n, eh) = 0;                   //set number of constraints acquired
        Vector3d constraint = Vector3d::Zero(3);    //a - constraint
        double bside = 0;                           //b - right side number

        std::set<MyMesh::FaceHandle> face_handles;
        std::set<MyMesh::VertexHandle> vertex_handles;
        double face_area;
        Vector3d normal = Vector3d::Zero(3);        //face normal coords storage
        MatrixXd v_coords = MatrixXd::Zero(3,3);    //matrix for vertex coordinates to compute determinant

        MatrixXd Hv = MatrixXd::Zero(3,3);          //Hessian for volume optimization
        MatrixXd Hb = MatrixXd::Zero(3,3);          //Hessian for boundary optimization
        MatrixXd Hs = MatrixXd::Zero(3,3);          //Hessian for triangle shape optimization
        Vector3d cv = Vector3d::Zero(3);            //for volume optimizaton
        Vector3d cb = Vector3d::Zero(3);            //for boundary optimization
        Vector3d cs = Vector3d::Zero(3);            //for triangle shape optimization
        double kv = 0;                              //constants in volume optimization
        double kb = 0;                              //constants in boundary optimization
        int i = 0;
        MatrixXd E1 = MatrixXd::Zero(3,3);          //e1 for every vertex
        MatrixXd E2 = MatrixXd::Zero(3,3);          //e2 for every vertex
        MatrixXd e1x = MatrixXd::Zero(3,3);         //e1x matrix for boundary optimization
        Vector3d e1 = Vector3d::Zero(3);            //summed E1
        Vector3d e2 = Vector3d::Zero(3);            //summed E2
        Vector3d e3 = Vector3d::Zero(3);            //cross product of e1 and e2
        Vector3d v0 = Vector3d::Zero(3);            //vertex coords storage
        Vector3d v1 = Vector3d::Zero(3);            //second vertex coords storage
        double fv = 0;                              //volume objective functions (the resulting error)
        double fb = 0;                              //area objective functions (the resulting error)

        Vector3d V = Vector3d::Zero(3);             //final vertex

    //---------------------------------------------------------------------------------------------------------------
        //Volume preservation (+volume optimization)
            //melo by jit pro vsechny, ale pokud by ta mesh byla az moc degenerovana, vyjde nula a v tom pripade ten constraint zahodime
            heh = mesh.halfedge_handle(eh, 0);          
            vh1 = mesh.to_vertex_handle(heh);           
            vh2 = mesh.from_vertex_handle(heh);
if(!(mesh.is_boundary(vh1) xor mesh.is_boundary(vh2))){
            //plug the faces we need into a set (every element is unique)
            for (auto vertex_handle : {vh1, vh2}) for (vf_it = mesh.vf_iter(vertex_handle); vf_it.is_valid(); ++vf_it) face_handles.insert(*vf_it);
            //calc constraint
            for (const auto& face_handle : face_handles) {
                if(mesh.is_valid_handle(face_handle)){
                    //compute face normal (should be magnitude 2x the area of the face) and the first constraint
                    face_normal = mesh.calc_face_normal(face_handle);
                    face_area = mesh.calc_face_area(face_handle);
                    face_normal = face_normal*2*face_area; 
                    for(int i = 0; i<3; i++) normal(i) = face_normal[i];
                    constraint += normal;
                    //get the determinant of the face and compute the first bside
                    for (i = 0, fv_it = mesh.fv_iter(face_handle); fv_it.is_valid(); ++fv_it, i++) {
                        p = mesh.point(*fv_it);
                        for (int j = 0; j<3; j++) v_coords(j, i) = p[j];
                    }
                    double determinant = v_coords.determinant();
                    bside += determinant;

                    //Volume optimization section
                    Hv += normal*normal.transpose();
                    cv -= determinant*normal.transpose();
                    kv += std::pow(determinant,2);

                }
            }
            if(is_alpha_compatible(eh, constraint)) add_constraint(eh, constraint, bside);


    //------------------------------------------------------------------------------------------------------------------------------
        //Boundary preservation
            //pouze pro hrany
        if(mesh.is_boundary(eh)){
            //std::cout<<"boundary found"<<std::endl;
            MyMesh::EdgeHandle* boundary_edges = new MyMesh::EdgeHandle[3];
            //get all needed handles
            heh = mesh.halfedge_handle(eh, 0); 
            eh1 = mesh.edge_handle(mesh.ccw_rotated_halfedge_handle(heh));
            eh2 = mesh.edge_handle(mesh.cw_rotated_halfedge_handle(heh));
            boundary_edges[0] = eh1; boundary_edges[1] = eh2; boundary_edges[2] = eh;

            //calculate E1 and E2
            for (int i = 0; i<3; i++) {
                eh = boundary_edges[i];
                heh = mesh.halfedge_handle(eh, 0);          
                vh1 = mesh.to_vertex_handle(heh);           
                vh2 = mesh.from_vertex_handle(heh);
                p0 = mesh.point(vh1);
                p1 = mesh.point(vh2);
                for (int j = 0; j<3; j++) {v0(j) = p0[j]; v1(j) = p1[j];}
                E1.col(i) = v1-v0;
                E2.col(i) = v1.cross(v0);
            }
            for (int i = 0; i<3; i++) {e1 += E1.col(i); e2 += E2.col(i);}   //sum up E1 and E2
            e3 = e1.cross(e2);
            delete[] boundary_edges;

            //equation 7
            constraint = e1.transpose()*e1*e3.transpose();
            bside = e3.transpose()*e3;
            if(is_alpha_compatible(eh, constraint)) add_constraint(eh, constraint, bside);
            //equation 8
            constraint = (e1.cross(e3)).transpose();
            bside = 0;
            if(is_alpha_compatible(eh, constraint)) add_constraint(eh, constraint, bside);
        }


    //-----------------------------------------------------------------------------------------------------------------------     
        //Volume optimization
            //da zbyvajici 2 constrainty pro ne-hrany
            //compute I and orthogonal Z (only needed if we don't have 3 constraints)
            if(mesh.property(n, eh) != 3) calc_remaining_constraints(eh, Hv, cv);

    //----------------------------------------------------------------------------------------------------------------------
        //Boundary optimization
        if(mesh.is_boundary(eh)){
            for (int i = 0; i<3; i++) {
                //create the (e x ) matrices
                e1x(0,1) = -E1.col(i)(2); e1x(0,2) = E1.col(i)(1); e1x(1,0) = E1.col(i)(2);
                e1x(1,2) = -E1.col(i)(0); e1x(2,0) = -E1.col(i)(1); e1x(2,1) = E1.col(i)(0);
                Hb += e1x.transpose()*e1x;
                temp1 = E1.col(i);
                temp2 = E2.col(i);
                cb += (temp1).cross(temp2);
                kb += (E2.col(i).transpose()*E2.col(i)).value();
            }
            //compute I and orthogonal Z (only needed if we don't have 3 constraints)
            if(mesh.property(n, eh) != 3) calc_remaining_constraints(eh, Hb, cb);
        }
        
    //----------------------------------------------------------------------------------------------------------------------
        //Apply triangle shape opt. if necessary
        if(mesh.property(n, eh) < 3){
            //insert needed vertices into a set
            for (auto vertex_handle : {vh1, vh2}) for (vv_it = mesh.vv_iter(vertex_handle); vv_it.is_valid(); ++vv_it) vertex_handles.insert(*vv_it);
            //and erase those, which are not needed
            vertex_handles.erase(vh1);
            vertex_handles.erase(vh2);
            //calculate the Hessian and cs
            for (auto vertex_handle : vertex_handles){
                    p = mesh.point(vertex_handle);
                    for (int j = 0; j<3; j++) {cs(j) -= p[j]; Hs(j,j)++;}                    
            }
            calc_remaining_constraints(eh, Hs, cs);
        }

        write_c(eh);

    //----------------------------------------------------------------------------------------------------------------------
        //Calculate edge collapse error
        if(mesh.property(n, eh) == 3){
            //get final vertex position
                mesh.property(v, eh).v = (mesh.property(a, eh).c.transpose()).colPivHouseholderQr().solve(mesh.property(b, eh).b);
                V = mesh.property(v, eh).v;
            //calculate volume and area optimization errors
                fv = (1/18.0)*((1/2.0)*(V.transpose()*(Hv*V)).value()+(cv.transpose()*V).value())+((1/36.0)*kv);
                fb = (1/2.0)*((1/2.0)*(V.transpose()*(Hb*V)).value()+(cb.transpose()*V).value())+((1/4.0)*kb);
            //calculate final error
                if(mesh.is_boundary(eh)){//use full formula
                    p0 = mesh.point(vh1);
                    p1 = mesh.point(vh2);
                    double length = (p1-p0).norm();
                    mesh.property(e, eh) = lambda*fv + (1-lambda)*std::pow(length,2)*fb;
                }
                else mesh.property(e, eh) = lambda*fv;
            //and push back to the vector of edges
                if(!init) eh_arr.push_back(eh);
        }

}
//else if(init) eh_arr.erase(std::find(eh_arr.begin(), eh_arr.end(), eh)); //if edge is a semiboundary edge (exactly one vertex is boundary vertex)


    }
}

//================================================================================================================================================
bool MeshWrap::is_alpha_compatible(MyMesh::EdgeHandle eh, Vector3d constraint){
    if(mesh.property(n, eh) == 0){
        if(constraint(0) == 0 and constraint(1) == 0 and constraint(2) == 0) return false;
        else return true;
    }
    else if(mesh.property(n, eh) == 1){
        if(std::pow((mesh.property(a, eh).c.col(0).transpose()*constraint), 2)
         < std::pow((mesh.property(a, eh).c.col(0).norm()*constraint.norm()*cos(alpha)),2)) return true;
        else return false;
        
    }
    else if(mesh.property(n, eh) == 2){
        temp1 = mesh.property(a, eh).c.col(0);
        temp2 = mesh.property(a, eh).c.col(1);
        temp3 = temp1.cross(temp2);
        if(std::pow((temp3.transpose()*constraint), 2)
         > std::pow((temp3.norm()*constraint.norm()*sin(alpha)),2)) return true;
        else return false;
    }
    return false;
}

void MeshWrap::add_constraint(MyMesh::EdgeHandle eh, Vector3d constraint, double right_side){
    mesh.property(a, eh).c.col(mesh.property(n, eh)) = constraint;
    mesh.property(b, eh).b(mesh.property(n, eh)) = right_side;
    mesh.property(n, eh)++;
}

void MeshWrap::calc_remaining_constraints(MyMesh::EdgeHandle eh, MatrixXd Hessian, Vector3d c){
    //create identity matrix of a size (3-n, 3) 
        MatrixXd I = MatrixXd::Zero(3-mesh.property(n, eh),3);
        for (int i = 2-mesh.property(n, eh), j = 2; i>=0; i--, j--) I(i, j) = 1;
    //create orthogonal matrix Z from known constraints
        MatrixXd Z = MatrixXd::Zero(3,3);
        for (int i = mesh.property(n, eh)-1; i>=0; i--) Z.col(i) = mesh.property(a, eh).c.col(i);   //plug in discovered constraints

        if(mesh.property(n, eh) == 0) Z = MatrixXd::Identity(3,3);
        else {
            if (mesh.property(n, eh) == 1) {Z(0,1) = Z(1,0); Z(1,1) = -Z(0,0);} 
            temp1 = Z.col(0);
            temp2 = Z.col(1);
            Z.col(2) = temp1.cross(temp2);     //add second orthogonal vector
        }
        
    //compute remaining constraints and b sides
        MatrixXd res = MatrixXd::Zero(3-mesh.property(n, eh),3);
        VectorXd resb = VectorXd::Zero(3-mesh.property(n, eh));
        auto temp = I*Z.inverse();
        res = temp*Hessian;
        resb = -(temp*c);
    //add constraints if possible
        for (int i = 2-mesh.property(n, eh); i>=0; i--) if(is_alpha_compatible(eh, res.row(i))) add_constraint(eh, res.row(i), resb(i));
           
}