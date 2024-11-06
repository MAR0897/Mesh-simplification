# Mesh simplification

## Obsah repozitáře (branches)

- Garland-Heckbert = samostatná aplikace nad OpenMesh knihovnou a implementace v OpenMesh systému
- GH-unconnected_vertices = pokus o implementaci spojování vrcholů v Garland-Heckbert simplifikaci
- Lindstrom-Turk = samostatná aplikace nad OpenMesh knihovnou
- LT-OpenMesh-system-impl = implementace Lindstrom-Turk simplifikace do decimačního systému OpenMeshe

### Spuštění LT simplifikace:
#### Lépe a rychleji to jede přes implementaci do OpenMesh systému v `LT-OpenMesh-system-impl`! 
cmake .\
make\
./main [input_file_path] [output_file_path] [number_of_vertices_to_be_simplified]

třeba:\
./main bunny.obj bunnyout.obj 200


