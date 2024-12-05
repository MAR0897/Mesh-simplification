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
- ```Decimater.cc```

### Spuštění
Stejný návod jako v LindTurk OpenMesh system implementation větvi

### Parametry
- lock boundary edges - true/false
- ideal vertex search mod
    - 0 (originální OpenMesh implementace, jediná změna je, že se přepočítává trochu více vrcholů (ale je to docela nevýznamné), ale to se dá v DecimaterT_impl.hh změnit)
    - 1 (počítání erroru pouze pro v0, v1, midpoint)
    - 2 (hledá se na přímce v0v1)
    - 3 (originální GH, hledá se v celém 3D prostoru)
- max error

### Command line spuštění
například:
```
    ./commandlineDecimater -i bunny.obj -o bunnyout.obj -M Q:true:3 -n 2000
```
zavolá originální GH decimaci a
```
    ./commandlineDecimater -i bunny.obj -o bunnyout.obj -M Q:true:0 -n 2000
```
decimuje mesh OpenMesh implementací

### Poznámky 
- GH přepočítává pouze pro sousední hrany result vertexu, aka pro ty, co se změní Q1+Q2
- numericky neoptimalizováno
- Všechny módy kromě originální OpenMesh implementace mají pro každý halfedge parametr, jestli už byl spočítán error opposite halfedge (tedy jestli má cenu to počítat znovu). Teoreticky by to tedy mělo být 2x rychlejší (pokud by tedy přepočet erroru bral násobně více času než všechno kolem toho, což zrovna u GH moc nenastává). U jiných algoritmů to ale může mít značnou úsporu.
- GH nepočítá s boundary edges (simplifikace bude probíhat špatně!!! - pomalu to bude ukousávat hranici meshe), takže se to musí pořešit. Mají tam své řešení, jinak se ty hrany musí locknout. V mé implementaci je možnost hrany locknout pomocí parametru.

### Přibližné výsledky
Collapses/s (počet odstraněných hran za sekundu):
- 0 = OpenMesh originální implementace = 50000
- 1 = výběr mezi v0, v1 a midpointem =   41000
- 2 = hledání na přímce v0v1 =           42000
- 3 = hledání v celém prostoru =         40000

Jde tedy vidět, že metody se blíží rychlosti originální OpenMesh implementace (i když jsou přesnější). 

Výhoda OpenMesh implementace je nižší počet hran, které se musí přepočítat (asi tak o počet sousedních hran odstraněného vrcholu). OpenMesh systém ale přepočítává vždy outgoing halfedge vrcholů, které jsou sousední s odstraňovaným vrcholem, takže ve výsledku je to ještě pomalejší, ale asi je to zase obecnější pro více druhů simplifikací. Kdyby byl ten systém postavý na hranách a ne halfedges, tak stačí přepočítat pouze sousední hrany zůstávajícího vrcholu. Ale dle papíru GH je to vlastně špatně to přepočítávání, protože zkolabovanému vrcholu se updatuje kvadrika, takže by se pro jeho sousední hrany měl přepočítat error, což se neděje (respektive děje, ale pouze na jednom ze 2 halfedgů, takže ten halfedge, pro který se error nepřepočítá, ho může mít nižší a bude to špatně) => upravil jsem kód, ať vybírá sousední vrcholy zůstávajícího vrcholu a ne toho, který se má odstranit.

Výhoda implementace dle GH je, že se nemusí počítat cost error pro oba halfedge (bo ideální pozice vrcholu, do kterého bychom chtěli hranu zkolabovat hledáme z celého prostoru, ne jenom jeden vrchol jako OpenMesh), ale stačí jeden (tedy kdyby systém byl postaven jenom na hranách, tak by to možná bylo i rychlejší). Ale kdyby se postavil systém pouze na hrany, tak ač by to bylo ještě rychlejší, zase moc nevidím, jaké jsou tradeoffy pak s tím sortováním cost hodnot. Náročnost sortování je nlog(n) a bavíme se o tom, jestli se bude sortovat vektor o velikosti počtu vrcholů nebo hran. U výpočtu erroru by se pouze jednalo o polovinu výpočtů při inicializaci erroru a pak třeba jen čtvrtinu výpočtů při rekalkulaci.  

**OpenMesh decimovací systém počítá cost (error) pro každý halfedge (= 2x počet hran), ale ve výsledku se do fronty přidá jenom nejnižší cost daného vrcholu (který se přidá mezi ty, které se mají přepočítat) a jednoho halfedge, takže celková fronta je ve skutečnosti o hodně menší, než kdybychom ji vytvořili jenom z hran.** Je tedy na zvážení, jestli by sortování fronty o velikosti počtu hran bylo nákladnější, než počítání o dost více hran, než které jsou nutné. GH totiž vyžaduje přepočítat pouze sousední hrany zkolabovaného výsledného vrcholu.

V mé implementaci u módů 1-3 je při inicializaci cost erroru přidělena bool hodnota všem hranám unikátním halfedgům, takže se skutečný výpočet probíhá pouze jednou u každé hrany. Systém ale stále cyklí přes halfedge a stále se může stát, že při přepočítávání se přepočítají i nějaké navíc, i třeba 3x tolik, než je potřeba.

Originální OpenMesh GH simplifikace je lepší na boundary hrany, protože je lépe zachovává. Jak již bylo řečeno, GH neřeší boundary hrany a proto je třeba je pořešit například locknutím. Když se nelocknou, tak hranice meshe se má tendenci smršťovat. OpenMesh implementace tak svou "nekorektností" vyřešila tento problém, protože u ní se to neděje.
