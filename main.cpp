#include "meshwrap.h"
                                                      
int main(int argv, const char **argc){//./main -i input_file_path -o output_file_path -n vertices_to_decimate
//./main -i objfiles/bunny.obj -o objfiles/bunnyout.obj -n 1 

    //Lyra arguments
    std::string input_path = argc[1];
    std::string output_path = argc[2];
    int number_of_vertices = std::atoi(argc[3]); 

    auto start = std::chrono::high_resolution_clock::now();

    MeshWrap m(input_path, output_path);
    std::cout << "Mesh successfully loaded into memory" << std::endl;
    std::cout<<"========================="<<std::endl;
    m.lock_boundary_edges();

    m.initialize();
    std::cout << "Mesh successfully initialized error on all edges" << std::endl;
    std::cout<<"========================="<<std::endl;

    m.simplify(number_of_vertices);
    std::cout<<"Mesh was successfully simplified"<<std::endl;
    std::cout<<"========================="<<std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    auto cas = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    std::cout<<"Took "<<cas<<" nanoseconds"<<std::endl;
        
    

return 0;}