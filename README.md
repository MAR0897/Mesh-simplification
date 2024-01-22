# Mesh-simplification

## Spuštění:
cmake .\
make\
./main -i [input_file_path] -o [output_file_path] -n [number_of_vertices_to_be_simplified]

třeba:\
./main -i objfiles/bunny.obj -o objfiles/bunnyout.obj -n 200

## Test 1 (1000/2500 vrcholu)
- Mesh successfully loaded into memory
- The code part #0 took 135 milliseconds to run.
- Mesh successfully initialized error on all edges 
- The code part #1 took 5564 milliseconds to run. 
- Mesh was successfully simplified
- The code part #2 took 34662 milliseconds to run.
- Cas funkce get_constraints_and_error: 852
- Cas funkce collapse_edge: 576
- Cas vybirani hrany na kolaps: 4928
- Cas vnejsiho cyklu rekalkulovani erroru u sousednich hran: 21933
- Cas pouze vnitrni funkce rekalkulovani erroru: 217

## Test 2 (2000/2500 vrcholu)
- Mesh successfully loaded into memory
- The code part #0 took 148 milliseconds to run.
- Mesh successfully initialized error on all edges
- The code part #1 took 5218 milliseconds to run.
- Mesh was successfully simplified
- The code part #2 took 64633 milliseconds to run.
- Cas funkce get_constraints_and_error: 2560
- Cas funkce collapse_edge: 745
- Cas vybirani hrany na kolaps: 7598
- Cas vnejsiho cyklu rekalkulovani erroru u sousednich hran: 47692
- Cas pouze vnitrni funkce rekalkulovani erroru: 2268

## Test 3 (2400/2500 vrcholu)
- Mesh successfully loaded into memory
- The code part #0 took 158 milliseconds to run.
- Mesh successfully initialized error on all edges
- The code part #1 took 5836 milliseconds to run.
- Mesh was successfully simplified
- The code part #2 took 75040 milliseconds to run.
- Cas funkce get_constraints_and_error: 3401
- Cas funkce collapse_edge: 739
- Cas vybirani hrany na kolaps: 7532
- Cas vnejsiho cyklu rekalkulovani erroru u sousednich hran: 57019
- Cas pouze vnitrni funkce rekalkulovani erroru: 2698


