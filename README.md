# LT simplification in OpenMesh decimater system

Zde se nacházejí nové nebo pozměněné zdrojové nebo hlavičkové soubory, které umožňují spustit Lindstrom-Turk simplifikaci meshe pomocí commandlineDecimator-u v OpenMesh systému

### Nové soubory

- `ModLindTurkT.hh` = deklarace funkcí pro LT simplifikaci (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `ModLindTurkT_impl.hh` = inicializace, výpočet erroru a preprocess_collapse (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)

### Pozměněné soubory

- `DecimaterT_impl.hh` = pozměněn způsob přepočítávání erroru (`OpenMeshRoot/src/OpenMesh/Tools/Decimater`)
- `decimater.cc` = hlavní aplikace pro decimaci v OpenMeshi, přidán mód pro LT simplifikaci (`OpenMeshRoot/src/OpenMesh/Apps/Decimating`)

### Spuštení

1. Stáhnout OpenMesh knihovnu a rozbalit
2. Vložit `ModLindTurkT.hh`, `ModLindTurkT_impl.hh` a `DecimaterT_impl.hh` do `OpenMeshRoot/src/OpenMesh/Tools/Decimater` a `decimater.cc` do `OpenMeshRoot/src/OpenMesh/Apps/Decimating`
3. V OpenMesh adresáři (OpenMeshRoot = OpenMesh-11.0.0 nebo podobně) zavolat:

    ```
    mkdir build
    cd build
    cmake ..
    make
    ```
5. Pro spuštění LT simplifikace stačí ve složce `build/Build/bin` zavolat:

    ```
    ./commandlineDecimater -i [input-file] -o [output-file] -M LT -n [n-of-vertices-to-decimate]
    ```
    Např.
    ```
    ./commandlineDecimater -i bunny.obj -o bunnyout.obj -M LT -n 30000
    ```
    K parametru LT se ještě dá přidat string `:[bool lock],[double lambda],[double alpha]` pro nastavení vlastních parametrů. 
    - Lock parametr zamkne boundary edge, aby hranice meshe zůstala stále stejná.
    - Lambda parametr je váha s jakou se počítá finální error jedné hrany.
    - Alfa parametr je úhel, do kterého se stěny meshe berou jako koplanární.
  
    Tedy například:

    ```
    ./commandlineDecimater -i bunny.obj -o bunnyout.obj -M LT:true,0.5,0.0174533 -n 30000
    ```
    Což jsou pro lambdu a alfu defaultní hodnoty

### Poznámky

- OpenMesh systém jede decimaci přes heap vrcholů (namísto hran). Pro každý vrchol se spočítá error pro každý (outgoing) halfedge, vybere se ten nejmenší a vrchol se pak zařadí do heapu podle velikosti erroru. Ve výsledku se pak pro každou hranu (edge) počítá error 2x (2x halfedge). Problém ale nastává při přepočítávání erroru po odstranění 1 hrany: přepočítávájí se pouze sousední vrcholy odstraněného vrcholu (více o odstranění hrany se dá najít v `CollapseInfoT.hh`). Což asi funguje pro simplifikaci, která vůbec neposunuje s vrcholy při odstranění hrany, ale pro LT už to moc nefunguje (simplifikovaná mesh vypadá podstatně hůře). Dá se to spravit přidáním sousedů vrcholu, který po odstranění hrany zůstane, ale stále se může stát, že u hrany, kde se spočítá pouze 1 ze 2 halfedgů (*), bude mít ten nepřepočítaný halfedge menší error, takže může dojít ke špatné simplikifaci. 

(*) tzn přepočítaný halfedge bude mít ten správný collapse error pro danou hranu, nepřepočítaný už bude mít špatný error

- Změna v souboru `DecimaterT_impl.hh` je ve funkci `decimate(...)`, kde se přepočítává error i u sousedních vrcholů vrcholu, který zůstává po odstranění hrany. Viz CTRL+F `ZMENA`.

- Garland-Heckbert simplifikace ani neposunuje vrcholy na jejich ideální pozice.

- Výpočet erroru by měl být správný, stačí se tudíž zamyslet nad jeho přepočtem (a pak na zrychlení kódu a zjednodušení maticových operací...).
