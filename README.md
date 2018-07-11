# AnalysisMetrologyITSU

1) Measure_Planarity.C

funzione principale: int Measure_Planarity()
variabili globali da settare (le altre si possono lasciare come sono):
- FileName: path + nome del file delle misure
- fMitutoyoFile: flag da settare true se il file in input è prodotto dalla Mitutoyo, false se è un file semplice con x,y,z in tre colonne
- nXcoord: numero di coordinate x della misura (a grandi linee, tipicamente 2 nel caso di metrologia di HS)
- xnom[nXcoord]: valori delle coordinate x nominali (x=+-15 mm nel caso di metrologia di HS)  
- fPlanarityWRTnominal: flag per scegliere se fare la planarità rispetto al piano nominale (=true) oppure il piano medio dei punti (=false)
- drawmodulespositions: flag per attivare o meno il disegno della posizione dei moduli (solo grafica)
- zmin, zmax: valori massimi e minimi di Zmisurato-Zpiano (da cambiare se i punti sono fuori scala o troppo compressi)


2) ComputeResidualsToNominalPositions.C

funzione principale: ComputeResidualsToNominalPositions()
variabili globali da settare:
- fMitutoyoFile: flag da settare true se il file in input è prodotto dalla Mitutoyo, false se è un file semplice con x,y,z in tre colonne
- FileName: path + nome del file delle misure
- FileNameNom: path + nome del file delle posizioni nominali ("OL_marker_nominal_positions_HS.dat")

3) ComparePadPositions.C

funzione principale: ComparePadPositions()
variabili globali da settare:
- infile1, infile2: path + nome dei file delle misure da confrontare. Nel caso si confronti una misura del SR212 con una del SR222, quella del 212 va in infile1
- fMitutoyoFile: flag da settare true se il file in input è prodotto dalla Mitutoyo, false se è un file semplice con x,y,z in tre colonne
- fChangeRS1: flag da settare true se il primo file contiene misure prese nel SR212, per fare la trasformazione del SR

4) FinalMarkerPositionExtrapolationAndQA.C

funzione principale: FinalMarkerPositionExtrapolationAndQA()
i file vengono caricati con una GUI. Per eseguire l'analisi cliccare su "DoMetrology", in basso a destra nella finestra. E' necessario avere il file delle posizioni nominali ("OL_marker_nominal_positions_HS.dat") nella stessa cartella da cui si esegue la macro

