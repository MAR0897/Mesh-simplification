#include "meshwrap.h"
                                                      
int main(int argv, const char **argc){//./main -i input_file_path -o output_file_path -n vertices_to_decimate
//./main -i objfiles/bunny.obj -o objfiles/bunnyout.obj -n 1 

    //Lyra arguments
    std::string input_path;
    std::string output_path;
    int number_of_vertices;
    bool show_help = false;
    auto cli
        = lyra::help(show_help)
        | lyra::opt(input_path, "input_path")
            ["-i"]["--input"]("Provide path for input .obj file")
        | lyra::opt(output_path, "output_path")
            ["-o"]["--output"]("Provide path for output .obj file")
        | lyra::opt(number_of_vertices, "number_of_vertices")
            ["-n"]["--number_of_vertices"]("How many vertices should the program decimate");
    auto result = cli.parse({argv, argc});
    if (show_help) {std::cout << cli << std::endl;return 0;}
    if (!result){std::cerr << "Error in command line: " << result.message() << std::endl;return 1;}
        
    auto start = std::chrono::high_resolution_clock::now();

    MeshWrap m(input_path, output_path);
    std::cout << "Mesh successfully loaded into memory" << std::endl; m.time(start);
    std::cout<<"========================="<<std::endl;
    m.initialize();
    std::cout << "Mesh successfully initialized error on all edges" << std::endl; m.time(start);
    std::cout<<"========================="<<std::endl;
    m.simplify(number_of_vertices);
    std::cout<<"Mesh was successfully simplified"<<std::endl; m.time(start);
    std::cout<<"========================="<<std::endl;

return 0;}