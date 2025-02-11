//Parameters for the coalescence simulation program: fastsimcoal2
2 samples to simulate:
//Population effective sizes (number of genes)
NPOP1 sfspool 0
NPOP2 sfspool 1
//Samples sizes and samples age Note that impossible to have inbreeding with odd numbers of samples, so we simplify the structure by regrouping pops
64 0 0.30
24 0 0.15
//Growth rates: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0 MIG21
MIG12 0
//Migration matrix 1
0 0
0 0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
3 historical event
TDIV 0 1 1 RESIZE3 0 1
TVAR3 1 1 1 RESIZE4 0 1
TVAR4 1 1 1 RESIZE5 0 1
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per gen recomb and mut rates
FREQ 1 0 1.7e-9 OUTEXP


