CPP_FILES = main.cpp

program: ${CPP_FILES}
	g++ ${CPP_FILES} -fsanitize=address -larmadillo -I/usr/local/include/OpenMesh -L/usr/local/lib -lOpenMeshCore -o main
