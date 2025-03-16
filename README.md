# Simplifikace implementované v OpenMesh systému

Zde se nacházejí nové nebo pozměněné zdrojové nebo hlavičkové soubory, které umožňují spustit Garland-Heckbert, Lindstrom-Turk a spektrální simplifikaci meshe pomocí commandlineDecimator-u v OpenMesh systému

### Nové/dost pozměněné soubory

- `ModQuadricT.hh` = deklarace funkcí pro GH simplifikaci (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `ModQuadricT_impl.hh` = inicializace a výpočet erroru GH simplifikace (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- - `ModLindTurkT.hh` = deklarace funkcí pro LT simplifikaci (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `ModLindTurkT_impl.hh` = inicializace a výpočet erroru LT simplifikace (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- - `ModSpectralT.hh` = deklarace funkcí pro spektrální simplifikaci (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `ModSpectralT_impl.hh` = inicializace a výpočet erroru spektrální simplifikace (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)

### Lehce pozměněné soubory

- `DecimaterT_impl.hh` = pozměněn způsob přepočítávání erroru (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `decimater.cc` = hlavní aplikace pro decimaci v OpenMeshi, přidán mód pro LT a spektrální simplifikaci (`OpenMeshRoot/src/OpenMesh/Apps/Decimating`)

### Spuštení

1. Stáhnout OpenMesh knihovnu a rozbalit
2. Vložit zdrojáky všech 3 simplifikací a `DecimaterT_impl.hh` do `OpenMeshRoot/src/OpenMesh/Tools/Decimater` a `decimater.cc` do `OpenMeshRoot/src/OpenMesh/Apps/Decimating`
3. V OpenMesh adresáři (OpenMeshRoot = OpenMesh-11.0.0 nebo podobně) zavolat:

    ```
    mkdir build
    cd build
    cmake ..
    make
    ```
5. Pro spuštění LT simplifikace stačí ve složce `build/Build/bin` zavolat:

    ```
    ./commandlineDecimater -i [input-file] -o [output-file] -M [simplification-mod-and-its-parameters] -n [n-of-vertices-to-decimate]
    ```
    Např.
    ```
    ./commandlineDecimater -i bunny.obj -o bunnyout.obj -M ML:0,0 -n 30000
    ```
    Simplifikační módy:
   - Garland-Heckber = Q
   - Lindstrom-Turk = ML (Memoryless)
   - spektrální = SP

   Parametry:
   1. algoritmus vybírání ideálního collapse vertexu - buď z celého 3D prostoru, nebo z přímky procházející danou hranou, nebo ze 3 bodů (v0,v1,midpoint)
      - pro GH je 0 originální OpenMesh implementace, 1 je vybírání ze 3 bodů, 2 je vybírání z přímky, 3 je přesný GH algoritmus
      - pro LT je 0 vybírání z celého prostoru a 1 ze přímky
   2. Lock parametr zamkne boundary edge, aby hranice meshe zůstala stále stejná. (=true nebo false)



### Eigen a Spectra knihovny CMake include změna

- commandlineDecimater
- DecimaterGui
- VDProgMesh/mkbalancedpm
