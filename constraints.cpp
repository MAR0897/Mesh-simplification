#include "meshwrap.h"

//Computes the constraits (result vertex coords) and error of the edge
void MeshWrap::get_constraints_and_error(MyMesh::EdgeHandle eh){

    //Take those edges, which are not locked
    if(!mesh.property(is_locked, eh)){
        heh = mesh.halfedge_handle(eh, 0);          
        vh1 = mesh.to_vertex_handle(heh);           
        vh2 = mesh.from_vertex_handle(heh);
        //Take those edges, which do not have exactly 1 boundary vertex
        if(!(mesh.is_boundary(vh1) xor mesh.is_boundary(vh2))){

            //set variables to zero
            i = 0;                      //iterator
            mesh.property(n, eh) = 0;   //set number of constraints acquired
            constraint.setZero();       //a - constraint
            bside = 0;                  //b - right side number
            face_handles.clear();       //set of face handles (to avoid repeating calculations)
            vertex_handles.clear();     //set of vertex handles (to avoid repeating calculations)
            normal.setZero();           //face normal coords storage
            v_coords.setZero();         //matrix for vertex coordinates to compute determinant
            Hv.setZero();               //Hessian for volume optimization
            Hb.setZero();               //Hessian for boundary optimization
            Hs.setZero();               //Hessian for triangle shape optimization
            cv.setZero();               //for volume optimizaton
            cb.setZero();               //for boundary optimization
            cs.setZero();               //for triangle shape optimization
            kv = 0;                     //constants in volume optimization
            kb = 0;                     //constants in boundary optimization
            E1.setZero();               //e1 for every vertex
            E2.setZero();               //e2 for every vertex
            e1x.setZero();              //e1x matrix for boundary optimization
            e1.setZero();               //summed E1
            e2.setZero();               //summed E2
            e3.setZero();               //cross product of e1 and e2
            v0.setZero();               //vertex coords storage
            v1.setZero();               //second vertex coords storage
            fv = 0;                     //volume objective functions (the resulting error)
            fb = 0;                     //area objective functions (the resulting error)
            V.setZero();                //final vertex

//---------------------------------------------------------------------------------------------------------------
        //Volume preservation (+volume optimization)    
            //plug the faces we need into a set (every element is unique)
            for (auto vertex_handle : {vh1, vh2}) for (vf_it = mesh.vf_iter(vertex_handle); vf_it.is_valid(); ++vf_it) face_handles.insert(*vf_it);
            //calc constraint
            for (const auto& face_handle : face_handles) {
                if(mesh.is_valid_handle(face_handle)){
                    //compute face normal (should be magnitude 2x the area of the face) and the first constraint
                    face_normal = mesh.calc_face_normal(face_handle);
                    face_area = mesh.calc_face_area(face_handle);
                    face_normal = face_normal*2*face_area;
                    Vector3d normal(face_normal[0], face_normal[1], face_normal[2]);
                    constraint += normal;
                    //get the determinant of the face and compute the first bside
                    for (i = 0, fv_it = mesh.fv_iter(face_handle); fv_it.is_valid(); ++fv_it, i++) {
                        p = mesh.point(*fv_it);
                        v_coords.col(i) = Vector3d(p[0], p[1], p[2]);
                    }
                    determinant = determinant3x3(v_coords);
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
            if(mesh.is_boundary(eh)){
                //get all needed handles
                MyMesh::EdgeHandle* boundary_edges = new MyMesh::EdgeHandle[3];
                heh = mesh.halfedge_handle(eh, 1); 
                eh1 = mesh.edge_handle(mesh.next_halfedge_handle(heh));
                eh2 = mesh.edge_handle(mesh.prev_halfedge_handle(heh));
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
                constraint = (e1.transpose()*e1)*e3.transpose();
                bside = -(e3.transpose()*e3).value();
                if(is_alpha_compatible(eh, constraint)) add_constraint(eh, constraint, bside);
                //equation 8
                constraint = (e1.cross(e3)).transpose();
                bside = 0;
                if(is_alpha_compatible(eh, constraint)) add_constraint(eh, constraint, bside);


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
            }

//-----------------------------------------------------------------------------------------------------------------------     
        //Volume optimization
            if(mesh.property(n, eh) != 3) calc_remaining_constraints(eh, Hv, cv);          
        
//----------------------------------------------------------------------------------------------------------------------
        //Boundary optimization
            if(mesh.is_boundary(eh) and mesh.property(n, eh) != 3) calc_remaining_constraints(eh, Hb, cb);

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

//----------------------------------------------------------------------------------------------------------------------
        //Calculate edge collapse error
            if(mesh.property(n, eh) == 3){
                //get final vertex position
                    mesh.property(v, eh).v = (mesh.property(a, eh).c.transpose()).colPivHouseholderQr().solve(mesh.property(b, eh).b);
                    V = mesh.property(v, eh).v;
                //calculate volume optimization error
                    fv = (V.transpose()*(Hv*V)).value() + 2*(cv.transpose()*V).value() + kv;
                //calculate final error
                    if(mesh.is_boundary(eh)){//use full formula
                        fb = 9*((V.transpose()*(Hb*V)).value() + 2*(cb.transpose()*V).value() + kb);
                        p0 = mesh.point(vh1);
                        p1 = mesh.point(vh2);
                        double length = (p1-p0).norm();
                        //std::cout<<Hv<<" "<<cb<<" "<<kb<<std::endl;
                        mesh.property(e, eh) = lambda*fv + (1-lambda)*std::pow(length,2)*fb;
                    }
                    else mesh.property(e, eh) = lambda*fv;
                //and push back to the vector of edges
                    if(!init) eh_arr.push_back(eh);
            }
        }
    }
}
//==================================================================================================================================================
//Checks if to-be-added constraint is alpha compatible
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

//Adds constraint to the system
void MeshWrap::add_constraint(MyMesh::EdgeHandle eh, Vector3d constraint, double right_side){
    mesh.property(a, eh).c.col(mesh.property(n, eh)) = constraint;
    mesh.property(b, eh).b(mesh.property(n, eh)) = right_side;
    mesh.property(n, eh)++;
}

//Calculates remaining constraints if there are less than 3
void MeshWrap::calc_remaining_constraints(MyMesh::EdgeHandle eh, MatrixXd Hessian, Vector3d c){
    //Create identity matrix of a size (3-n, 3)
        MatrixXd I = MatrixXd::Zero(3-mesh.property(n, eh),3);
        for (int i = 2-mesh.property(n, eh), j = 2; i>=0; i--, j--) I(i, j) = 1;
    //Create orthogonal matrix Z
        MatrixXd Z = MatrixXd::Zero(3,3);
        for (int i = mesh.property(n, eh)-1; i>=0; i--) Z.col(i) = mesh.property(a, eh).c.col(i);
        if(mesh.property(n, eh) == 0) Z = MatrixXd::Identity(3,3);
        else {
            if (mesh.property(n, eh) == 1) {Z(0,1) = Z(1,0); Z(1,1) = -Z(0,0);}     //Add first orthogonal vector
            temp1 = Z.col(0);
            temp2 = Z.col(1);
            Z.col(2) = temp1.cross(temp2);                                          //Add second orthogonal vector
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