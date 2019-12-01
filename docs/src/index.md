# InformMe.jl Documentation



```@contents
```

## Functions


```@docs
fastaToCpG(FASTAfilename;
                    outdir="./",
                    wsize=1000)
```

```@docs
convertBAMtoBits(bamFilename,phenoName,chr_num;
                          reference_path="./genome/",
                          bamfile_path="./indexedBAMfiles/",
                          matrices_path="./matrices/",
                          estimation_path="./estimation/",
                          scratch_path="./scratch/",
                          outdir="./output/",
                          pairedEnds=true,
                          numBasesToTrim=0,
                          minCpGsReqToModel=10,
                          regionSize=3000,
                          boundaryConditions=false,
                          MSIflag=false,
                          ESIflag=false,
                          MCflag=false,
                          subRegionSize=150)
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
