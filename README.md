# GH Mesh-simplification

## Custom implementace
soubory:
- ```main.cpp```
- ```CMakelists.txt```
- ```lyra.hpp```

### Custom implementace - Spuštění:
cmake .\
make\
./main [input_file_path] [output_file_path] [number_of_vertices_to_be_simplified] [decimator_mode]

### Custom implementace - Decimator mode
1 = OpenMesh inbuilt quadric decimator\
2 = moje implementace (to samé akorát s locked edges a o hodně pomalejší)

## Implementace v OpenMesh systému
soubory:
- DecimaterT_impl.hh
- ModQuadricT.hh
- ModQuadricT_impl.hh

### Spuštění
Stejný návod jako v LindTurk OpenMesh system implementation větvi

### Parametry
- lock boundary edges - true/false
- ideal vertex search mod - 0 (originální OpenMesh implementace, jediná změna je, že se přepočítává trochu více vrcholů, ale to se dá v DecimaterT_impl.hh změnit), 1 (počítání erroru pouze pro v0, v1, midpoint), 2 (hledá se na přímce v0v1), 3 (originální GH, hledá se v celém 3D prostoru)
- max error

### Poznámky 
- GH přepočítává pouze pro sousední hrany result vertexu, aka pro ty, co se změní Q1+Q2
- numericky neoptimalizováno

