# InformMe.jl Documentation

```@contents
```

## RUNNING informME

Run the informME software using the following steps in the indicated order (see below for the detailed documentation of each function).

### 1. REFERENCE GENOME ANALYSIS:
	
`fastaToCpg(FASTAfilename)`

This step analyzes the reference genome FASTA\_FILE (in FASTA format) and produces a julia JLD2 file CpGlocationChr#.jld2 for each chromosome, which is stored by default in REFGENEDIR, and contains the following information:

* location of CpG sites 

* CpG density for each CpG site 

* distance between neighboring CpG sites

* location of the last CpG site in the chromosome

* length of chromosome (in base pairs)

NOTE1: This step only needs to be completed one time for a given reference genome. Start analyzing samples at step 2 if you have previously completed step 1 for your sample's reference genome.

NOTE2: At this time the statistical model of informME has been designed to work only with autosomes, and so the informME software is not recommended for mitochondrial chromosomes, lambda spike-ins, partial contigs, sex chromosomes, et cetera. 

### 2-5. COMBINED SINGLE SAMPLE ANALYSIS:

`convertBAMtoBits(bamFilenames,phenoName)`

Steps 2 through 5 below can be invoked by the single command above, where bamFilenames is a list of bam files to be modeled and assigned the phenotype name given by the argument phenoName which should be unique to this sample. 


### 2. METHYLATION DATA MATRIX GENERATION: 
	
`matrixFromBam(BAM_FILE,CHR_NUM)`

This step takes the BAM file BAM\_FILE as input and generates the methylation data matrix for chromosome number CHR\_NUM. By default, the file BAM\_FILE and its associated index file (with extension .bai) is expected to be in BAMDIR, and the output file produced by this step is stored in a subdirectory in INTERDIR named after the chromosome number CHR\_NUM. The output file preserves the prefix from the file BAM\_FILE and the suffix '\_matrices.mat' is appended to it (e.g. if BAM\_FILE is normal\_sample.bam and CHR\_NUM is 10, then the output file is saved as INTERDIR/chr10/normal\_sample\_matrices.mat). The file produced contains the following information for each genomic region, which is subsequently used for model estimation:

* data matrix with -1,0,1 values for methylation status

* CpG locations broken down by region

NOTE1: See reference [1], "Online Methods: Quality control and alignment" for our suggested preprocessing steps when generating a sorted, indexed, deduplicated BAM file to input to informME.  

### 3. MODEL ESTIMATION:

`estParamsForChr(MAT_FILES,PHENO,matrices_path,reference_path,CHR_NUM)`

informME learns the parameters of the Ising probability distribution by combining the methylation data matrices provided through the argument MAT\_FILES (comma-separated list) for chromosome number CHR\_NUM. By default, the MAT\_FILES are expected to be in a subdirectory named after CHR\_NUM in "./results/". The output generated during this phase is also stored in a subdirectory in results named after chromosome number CHR\_NUM. The output file has as prefix PHENO and the suffix '\_fit.jld2' appended to it (e.g. if 'normal' is the PHENO, and CHR\_NUM is 10, then the output is stored as ./results/chr10/normal\_fit.jld2). The file produced contains the following information:

* CpG distances

* CpG densities

* estimated alpha, beta, and gamma parameters of the Ising model

* initial and transition probabilities of the inhomogeneous Markov chain representation of the Ising model

* marginal probabilities at each CpG site

* the log partition function of the estimated Ising model

### 4. MODEL ANALYSIS:

`methAnalysisForChr(prefix,chr_num,reference_path,estimation_path)`

This step consists of analyzing the model learned by computing a number of statistical summaries of the methylation state, including probability distributions of methylation levels, mean methylation levels, and normalized methylation entropies, as well as mean and entropy based classifications. This step also computes entropic sensitivity indices, methylation sensitivity indices, as well information-theoretic quantities associated with methylation channels, such as turnover ratios, channel capacities, and relative dissipated energies. The output generated during this phase is stored in the same directory as the output generated during the first phase, using the same prefix as before. However, the suffix is now '\_analysis.mat' (e.g. following the previous example, the output file of this phase is stored as ./results/chr10/normal\_analysis.mat). This file contains the following information:

* the locations of the CpG sites within the genomic region

* numbers of CpG sites within the analysis subregions 

* which analysis subregions are modeled and which are not

* estimated parameters of Ising model in genomic region

* methylation level probabilities in modeled subregions

* coarse methylation level probabilities

* mean methylation levels

* normalized methylation entropies

* entropic sensitivity indices 

* methylation sensitivity indices

* turnover ratios

* channel capacities

* relative dissipated energies


### 5. GENERATE BED FILES FOR SINGLE ANALYSIS:

`makeBedsForMethAnalysis(PHENO,analysis_path,reference_path)`

This function makes BED files from the methylation analysis results obtained after running methAnalysisForChr for a given phenotype PHENO. By default, the input file (analysis file) is expected to be located in ./results/chr#/PHENO\_analysis.mat. In addition, the output files are stored in "./singleMethAnalysisToBed_out/" and have the following names and content:

* MML-PHENO.bed: mean methylation levels

* NME-PHENO.bed: normalized methylation entropy
    
* METH-PHENO.bed: methylation-based classification (non-variable)
    
* VAR-PHENO.bed: methylation-based classification (variable)
    
* ENTR-PHENO.bed: entropy-based classification
    
* ESI-PHENO.bed (if ESIflag passed): entropic sensitivity indices

* MSI-PHENO.bed (if MSIflag passed): methylation sensitivity indices
    
* TURN-PHENO.bed (if MCflag passed): turnover ratios
    
* CAP-PHENO.bed (if MCflag passed): channel capacities
    
* RDE-PHENO.bed (if MCflag passed): relative dissipated energies  


### 6. GENERATE BED FILES FOR DIFFERENTIAL ANALYSIS:

`makeBedsForDiffMethAnalysis(PHENO1,PHENO2,analysis_path_1,
                             analysis_path_2,reference_path)`


This function makes BED files for the differential methylation analysis results obtained after running methAnalysisForChr for two given phenotypes PHENO1 and PHENO2. The input files (both analysis files) are expected to be located in analysis_path_1/chr#/PHENO1\_analysis.jld2 and analysis_path_2/chr#/PHENO2\_analysis.jld2 respectively. In addition, the output files are stored in "./makeBedsForDiffMethAnalysis_out/" and have the following names and content:

* dMML-PHENO1-VS-PHENO2.bed: differences in mean methylation levels
       
* DMU-PHENO1-VS-PHENO2.bed: differential mean-based classification
       
* dNME-PHENO1-VS-PHENO2.bed: differences in normalized methylation entropies
       
* DEU-PHENO1-VS-PHENO2.bed: differential entropy-based classification
       
* JSD-PHENO1-VS-PHENO2.bed: Jensen-Shannon distances
       
* dESI-PHENO1-VS-PHENO2.bed (if --ESI flag passed): differences in entropic sensitivity indices

* dMSI-PHENO1-VS-PHENO2.bed (if --MSI flag passed): differences in methylation sensitivity indices
       
* dCAP-PHENO1-VS-PHENO2.bed (if --MC flag passed): differences in channel capacities
       
* dRDE-PHENO1-VS-PHENO2.bed (if --MC flag passed): differences in relative dissipated energies


### 7. POST PROCESSING:

See the [original informME package](https://github.com/GarrettJenkinson/informME) for post processing scripts in R for gene/region rankings and DMR finding. 


## Functions


```@docs
fastaToCpG(FASTAfilename;
                    outdir="./",
                    wsize=1000)
```

```@docs
convertBAMtoBits(bamFilenames,phenoName;
                          reference_path="./genome/",
                          bamfile_path="./indexedBAMfiles/",
                          matrices_path="./matrices/",
                          estimation_path="./estimation/",
                          outdir="./output/",
                          pairedEnds=true,
                          numBasesToTrim=0,
                          minCpGsReqToModel=10,
                          regionSize=3000,
                          boundaryConditions=false,
                          MSIflag=false,
                          ESIflag=false,
                          MCflag=false,
                          subRegionSize=150,
                          chr_nums=1:22,
                          numProcessors=nprocs())
```


```@docs
matrixFromBam(bamFilename,chr_num;
                       reference_path="./genome/",
                       bamfile_path="./indexedBAMfiles/",
                       matrices_path="./matrices/",
                       pairedEnds=true,
                       numBasesToTrim=0,
                       regionSize=3000,
                       minCpGsReqToModel=10)
```


```@docs
estParamsForChr(mat_files,phenoName,matrices_path,reference_path,chr_num;
                         estimation_path="./results/",
                         regionSize=3000,
                         boundaryConditions=false)
```

```@docs
 methAnalysisForChr(phenoName,chr_num,reference_path,estimation_path;
                        outdir="./results/",
                        MSIflag=false,
                        ESIflag=false,
                        MCflag=false,
                        regionSize=3000,
                        subRegionSize=150)
```

```@docs
makeBedsForMethAnalysis(phenoName,analysis_path,reference_path;
                              outdir="./singleMethAnalysisToBed_out/",
                              chrs=[string("chr",i) for i=1:22],
                              #minChrNum=1,
                              #maxChrNum=22,
                              MSIflag=false,
                              ESIflag=false,
                              MCflag=false,
                              thresh=0.4,
                              regionSize=3000,
                              subregionSize=150)
```

```@docs
diffMethAnalysisToBed(phenoName_1,phenoName_2,analysis_path_1,analysis_path_2,reference_path,
                                  outdir="../makeBedsForDiffMethAnalysis_out/",
                                  chrs=[string("chr",i) for i=1:22],
                                  #  minChrNum=1,
                                  #  maxChrNum=22,
                                  MSIflag=false,
                                  ESIflag=false,
                                  MCflag=false,
                                  regionSize=3000,
                                  subregionSize=150,
                                  minNumCpG=2,
                                  thresh=0.55,
                                  threshDMU=[-1,-0.55,-0.1,0.1,0.55,1],
                                  threshDEU=[-1,-.5,-.3,-.05,.05,.3,.5,1])
```

## Index

```@index
```
