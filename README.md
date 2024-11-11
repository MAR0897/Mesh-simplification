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
- ```DecimaterT_impl.hh```
- ```ModQuadricT.hh```
- ```ModQuadricT_impl.hh```

### Spuštění
Stejný návod jako v LindTurk OpenMesh system implementation větvi

### Parametry
- lock boundary edges - true/false
- ideal vertex search mod
    - 0 (originální OpenMesh implementace, jediná změna je, že se přepočítává trochu více vrcholů, ale to se dá v DecimaterT_impl.hh změnit)
    - 1 (počítání erroru pouze pro v0, v1, midpoint)
    - 2 (hledá se na přímce v0v1)
    - 3 (originální GH, hledá se v celém 3D prostoru)
- max error

### Poznámky 
- GH přepočítává pouze pro sousední hrany result vertexu, aka pro ty, co se změní Q1+Q2
- numericky neoptimalizováno
- Všechny módy kromě originální OpenMesh implementace mají pro každý halfedge parametr, jestli už byl spočítán error opposite halfedge (tedy jestli má cenu to počítat znovu. Teoreticky by to tedy mělo být 2x rychlejsí (pokud by tedy přepočet erroru bral násobně více času než všechno kolem toho, což zrovna u GH moc nenastává). U jiných algoritmů to ale může mít značnou úsporu.
- GH nepočítá s boundary edges (simplifikace bude probíhat špatně!!!), takže se to musí pořešit. Mají tam své řešení, jinak se ty hrany musí locknout. V mé implementaci je možnost hrany locknout pomocí parametru.

### Výsledky
Collapses/s:
- 0 = 50000
- 1 = 41000
- 2 = 42000
- 3 = 40000

Jde tedy vidět, že i když jednotlivé módy nejsou numericky optimalizované, tak stejně se blíží originální OpenMesh implementaci. Výhoda OpenMesh implementace je nízký počet hran, které se musí přepočítat (mělo by se rovnat počtu sousedních hran vrcholu, který odstraňujeme (OpenMesh systém ale přepočítává vždy outgoing halfedge vrcholů, které jsou sousední s odstraňovaným vrcholem, takže ve výsledku je to ještě pomalejší, ale asi je to zase obecnější pro více druhů simplifikací)). Výhoda implementace dle GH je, že se nemusí počítat cost error pro oba halfedge, ale stačí jeden (tedy kdyby systém byl postaven jenom na hranách, tak by to možná bylo i rychlejší). GH implementace ale musí přepočítávat i sousední hrany vrcholu, který zůstává, takže přece jenom o něco více hran. Kdyby se ale postavil systém pouze na hrany, tak by to bylo ještě rychlejší, ale zase moc nevidím, jaké jsou tradeoffy pak s tím sortováním cost hodnot. **OpenMesh decimovací systém totiž počítá cost (error) pro každý halfedge (= 2x počet hran), ale ve výsledku se do fronty přidá jenom nejnižší cost daného vrcholu (který se přidá mezi ty, které se mají přepočítat) a jednoho halfedge, takže celková fronta je ve skutečnosti o hodně menší, než kdybychom ji vytvořili jenom z hran.** Je tedy na zvážení, jestli by sortování fronty o velikosti počtu hran bylo nákladnější, než počítání o dost více hran, než které jsou nutné. GH totiž vyžaduje přepočítat pouze sousední hrany zkolabovaného výsledného vrcholu.

V mé implementaci u módů 1-3 je při inicializaci cost erroru přidělena bool hodnota všem hranám unikátním halfedgům, takže se skutečný výpočet probíhá pouze jednou u každé hrany. Systém ale stále cyklí přes halfedge.
