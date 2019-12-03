# informME: An information-theoretic pipeline for WGBS data
# Copyright (C) 2018, Garrett Jenkinson (jenkinson@jhu.edu),
# and Jordi Abante (jabante1@jhu.edu)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
# or see <http://www.gnu.org/licenses/>.
#
@doc raw"""
 This function makes BED files for the methylation analysis results
 obtained by means of MethAnalysisForChr.m for a single phenotype.

## USAGE (default):

 `makeBedsForMethAnalysis(prefix,analysis_path,reference_path)`

## USAGE (optional):

 Example of optional usage with additional input parameters.

`makeBedsForMethAnalysis(prefix,analysis_path,reference_path,
                         outdir="/path/to/output")`

## MANDATORY INPUTS:

 `prefix`

               A string that specifies the name of the phenotype.

 `analysis_path`

               A string that specifies the path of the directory in which
               the model was constructed.

 `reference_path`

               A string that specifies the path to the directory that
               contains the results of analysis of the reference genome
               performed by FastaToCpG.m as well as the results of
               methylation calling performed by matrixFromBam.jl.

## OPTIONAL INPUTS:

`chrs`

              A vector of strings for the chromosomes to output to the
              final bed files. Default value: `[string("chr",i) for i=1:22]`

 `outdir`

               A string that specifies the path of the directory in which
               the output BED files are written.
               Default value: "./"


 `MSIflag`

               Flag that determines whether this function performs
               computation of the methylation sensitivity index (MSI).
               false: no MSI computation.
               true: allow MSI computation.
               Default value: false

 `ESIflag`

               Flag that determines whether this function performs
               computation of the entropic sensitivity index (ESI).
               false: no ESI computation.
               true: allow ESI computation.
               Default value: false

 `MCflag`

               Flag that determines whether this function performs
               computation of turnover ratios, CpG entropies, capacities,
               and relative dissipated energies of methylation
               channels (MCs).
               false: no MC computations.
               true: allow MC computations.
               Default value: false

 `thresh`

               A scalar used as a threshold in methylation-based
               classification.
               Default value: 0.4

 `regionSize`

               The size of the genomic regions used for parameter
               estimation (in number of base pairs).
               Default value: 3000

 `subregionSize`

               The size of the subregions of a genomic region used
               for methylation analysis (in number of base pairs).
               The ratio regionSize/subregionSize must be an integer.
               Default value: 150

 The default values of thresh, regionSize, and subregionSize should only
 be changed by an expert with a detailed understanding of the code and
 the methods used.
"""
function makeBedsForMethAnalysis(phenoName,analysis_path,reference_path;
                              outdir="./singleMethAnalysisToBed_out/",
                              chrs=[string("chr",i) for i=1:22],
                              MSIflag=false,
                              ESIflag=false,
                              MCflag=false,
                              thresh=0.4,
                              regionSize=3000,
                              subregionSize=150)


  # Manual checks/corrections of inputs
  if analysis_path[end:end]!="/"
    analysis_path=string(analysis_path, "/")
  end
  if reference_path[end:end]!="/"
    reference_path=string(reference_path, "/")
  end
  if outdir[end:end]!="/"
    outdir=string(outdir, "/")
  end

  #
  ###########################################################################
  # Hard coded options
  ###########################################################################
  # MML
  hardCodedOptionsMML = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=MML-"
  # NME
  hardCodedOptionsNME = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=NME-"
  # METH
  hardCodedOptionsMETH = "track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-2.0:2.0 name=METH-"
  # VAR
  hardCodedOptionsVAR = "track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=0.0:3.0 name=VAR-"
  # ENTR
  hardCodedOptionsENTR = "track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-2.0:2.0 name=ENTR-"

  # MSI
  MSIflag && (hardCodedOptionsMSI = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=MSI-")

  # ESI
  ESIflag && (hardCodedOptionsESI = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=ESI-")

  # TURN
  MCflag && (hardCodedOptionsTURN = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=TURN-")

  # CAP
  MCflag && (hardCodedOptionsCAP = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=CAP-")

  # RDE
  MCflag && (hardCodedOptionsRDE = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on alwaysZero=on name=RDE-")

  # Open .bed files to write output

  #check that bedfile folder exists, create it if not
  if !isdir(outdir)
    mkdir(outdir)
  end

  #if minChrNum > maxChrNum
  #  println("error: minChrNum must be <= maxChrNum")
  #  return
  #end

  # MML
  bedFileNameMML = string(outdir, "MML-", phenoName, ".bed")
  bedFileIDMML   = open(bedFileNameMML,"w")
  # NME
  bedFileNameNME = string(outdir, "NME-", phenoName, ".bed")
  bedFileIDNME   = open(bedFileNameNME,"w")
  # METH
  bedFileNameMETH = string(outdir, "METH-", phenoName, ".bed")
  bedFileIDMETH   = open(bedFileNameMETH,"w")
  # VAR
  bedFileNameVAR = string(outdir, "VAR-", phenoName, ".bed")
  bedFileIDVAR   = open(bedFileNameVAR,"w")
  # ENTR
  bedFileNameENTR = string(outdir, "ENTR-", phenoName, ".bed")
  bedFileIDENTR   = open(bedFileNameENTR,"w")

  if MSIflag
    # MSI
    bedFileNameMSI = string(outdir, "MSI-", phenoName, ".bed")
    bedFileIDMSI   = open(bedFileNameMSI,"w")
  end

  if ESIflag
    # ESI
    bedFileNameESI = string(outdir, "ESI-", phenoName, ".bed")
    bedFileIDESI   = open(bedFileNameESI,"w")
  end

  if MCflag
    # TURN
    bedFileNameTURN = string(outdir, "TURN-", phenoName, ".bed")
    bedFileIDTURN   = open(bedFileNameTURN,"w")
    # CAP
    bedFileNameCAP = string(outdir, "CAP-", phenoName, ".bed")
    bedFileIDCAP   = open(bedFileNameCAP,"w")
    # RDE
    bedFileNameRDE = string(outdir, "RDE-", phenoName, ".bed")
    bedFileIDRDE   = open(bedFileNameRDE,"w")
  end

  # write header info
  # MML
  @printf(bedFileIDMML,"%s%s\n",hardCodedOptionsMML,phenoName)
  # NME
  @printf(bedFileIDNME,"%s%s\n",hardCodedOptionsNME,phenoName)
  # METH
  @printf(bedFileIDMETH,"%s%s\n",hardCodedOptionsMETH,phenoName)
  # VAR
  @printf(bedFileIDVAR,"%s%s\n",hardCodedOptionsVAR,phenoName)
  # ENTR
  @printf(bedFileIDENTR,"%s%s\n",hardCodedOptionsENTR,phenoName)

  if MSIflag
    # MSI
    @printf(bedFileIDMSI,"%s%s\n",hardCodedOptionsMSI,phenoName)
  end

  if ESIflag
    # ESI
    @printf(bedFileIDESI,"%s%s\n",hardCodedOptionsESI,phenoName)
  end

  if MCflag
    # TURN
    @printf(bedFileIDTURN,"%s%s\n",hardCodedOptionsTURN,phenoName)
    # CAP
    @printf(bedFileIDCAP,"%s%s\n",hardCodedOptionsCAP,phenoName)
    # RDE
    @printf(bedFileIDRDE,"%s%s\n",hardCodedOptionsRDE,phenoName)
  end

  # brings data type into scope to help compiler
  regStruct::AnalyRegStruct = AnalyRegStruct()

  ###########################################################################
  # Loop over chromosomes
  ###########################################################################
  for chr_str in chrs
    #try
      if (length(chr_str)>3) && (chr_str[1:3]=="chr")
        chr_num_str = chr_str[4:end]
      else
        chr_num_str = chr_str
        chr_str = string("chr",chr_num_str)
      end

      #
      # Find last CpG site on chromosome (to determine how many regions to loop through)
      #
      # Load CpG site data on chromosome (to determine how many regions to loop through, etc.)
      CpGdata = string(reference_path, "CpGlocationChr", chr_num_str, ".jld2")
      finalCpGloc = FileIO.load(CpGdata,"finalCpGloc")
      CpGlocation = FileIO.load(CpGdata,"CpGlocation")

#       CpGdata = string(reference_path, "CpGlocationChr", chr_num_str, ".mat")
#       matfile=matopen(CpGdata)
#       finalCpGloc = read(matfile,"finalCpGloc")
#       CpGlocation = read(matfile,"CpGlocation")
#       close(matfile)

      #
      # find all start base-pairs for regions to be modeled on this chromsome
      #
      allStartBPs = 1:regionSize:finalCpGloc  # only need to go to the last CpG site, not last BP

      #
      # Load relevant results files
      #
      analysFile = string(analysis_path, chr_str, "/", phenoName, "_analysis.jld2")
      mapObjAnaly = FileIO.load(analysFile,"mapObjAnaly")

      #
      ###########################################################################
      # Loop through all regions on chromosome
      ###########################################################################
      #

      for regionNum = 1:length(allStartBPs)
        #try

          # find path of this region
          startBP = allStartBPs[regionNum]
          endBP   = startBP+regionSize-1
          locationPathName = string(chr_str, "/bp", startBP, "-", endBP)

          if haskey(mapObjAnaly, locationPathName)
            #
            # load structs for region
            #
            regStruct  = mapObjAnaly[locationPathName]
            cLprobs    = regStruct.cLprobs
            isModeled  = regStruct.isModeled
            Ncg        = regStruct.Ncg
            MML        = regStruct.MML
            NME        = regStruct.NME

            if MSIflag
              MSI = regStruct.MSI
              if isempty(MSI)
                  close(bedFileIDMML)
                  close(bedFileIDMETH)
                  close(bedFileIDVAR)
                  close(bedFileIDENTR)
                  close(bedFileIDNME)
                  close(bedFileIDMML)
                  close(bedFileIDMSI)
                  ESIflag && close(bedFileIDESI)
                  MCflag && close(bedFileIDTURN)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: MSIflag = true, but the previous methylation analysis step was run with MSIflag = false.\nRerun this julia function with MSIflag = false or rerun the previous methylation analysis step with MSIflag = true after deleting $analysFile.\n")
                  return 1
              end
            end

            if ESIflag
              ESI = regStruct.ESI
              if isempty(ESI)
                  close(bedFileIDMML)
                  close(bedFileIDMETH)
                  close(bedFileIDVAR)
                  close(bedFileIDENTR)
                  close(bedFileIDNME)
                  close(bedFileIDMML)
                  close(bedFileIDESI)
                  MSIflag && close(bedFileIDMSI)
                  MCflag && close(bedFileIDTURN)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: ESIflag = true, but the previous methylation analysis step was run with ESIflag = false.\nRerun this julia function with ESIflag = false or rerun the previous methylation analysis step with ESIflag = true after deleting $analysFile.\n")
                  return 1
              end
            end

            if MCflag
              TURN = regStruct.TURN
              CAP  = regStruct.CAP
              RDE  = regStruct.RDE
              if isempty(CAP) || isempty(RDE) || isempty(TURN)
                close(bedFileIDMML)
                  close(bedFileIDMETH)
                  close(bedFileIDVAR)
                  close(bedFileIDENTR)
                  close(bedFileIDNME)
                  close(bedFileIDMML)
                  MSIflag && close(bedFileIDMSI)
                  ESIflag && close(bedFileIDESI)
                  MCflag && close(bedFileIDTURN)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: MCflag = true, but the previous methylation analysis step was run with MCflag = false.\nRerun this julia function with MCflag = false or rerun the previous methylation analysis step with MCflag = true after deleting $analysFile.\n")
                  return 1
              end
            end

            #
            # Compute regional quantities
            #

            lower_index,upper_index = findSortedIndices(CpGlocation,startBP,endBP)
            CpGlocInReg = CpGlocation[lower_index:upper_index]
            numCpGinReg = length(CpGlocInReg)

            subRegCount = 0
            for offset = 0:subregionSize:(regionSize-1)
              subRegCount = subRegCount + 1
              subStartBP  = startBP + offset
              subEndBP    = min(subStartBP + subregionSize - 1,finalCpGloc)

              numCpGsInSubReg = Ncg[subRegCount]

              if numCpGsInSubReg>0 && isModeled[subRegCount] # region modeled

                lower_index,upper_index = findSortedIndices(CpGlocInReg,subStartBP,subEndBP)
                if numCpGsInSubReg==1 && (lower_index==1||upper_index==numCpGinReg) # only 1 cpg site and its a boundary condition
                  continue  #move onto next subregion
                end

                cLprobsRegion = cLprobs[subRegCount,:]
                MMLRegion     = MML[subRegCount]
                NMERegion     = NME[subRegCount]
                MSIflag && (MSIRegion = MSI[subRegCount])
                ESIflag && (ESIRegion = ESI[subRegCount])
                MCflag && (TURNRegion = TURN[subRegCount])
                MCflag && (CAPRegion  = CAP[subRegCount])
                MCflag && (RDEregion  = RDE[subRegCount])


                # Print out  MML,ESI,TURN,CAP,RDE

                MMLRegion<Inf &&  @printf(bedFileIDMML,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),MMLRegion)
                MSIflag && (MSIRegion<Inf) && @printf(bedFileIDMSI,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),MSIRegion)
                ESIflag && (ESIRegion<Inf) && @printf(bedFileIDESI,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),ESIRegion)
                MCflag && (TURNRegion<Inf) && @printf(bedFileIDTURN,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),CAPRegion)
                MCflag && (CAPRegion<Inf) && @printf(bedFileIDCAP,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),CAPRegion)
                MCflag && (RDEregion<Inf) && @printf(bedFileIDRDE,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),RDEregion)


                # Do NME classification

                if NMERegion<Inf
                  @printf(bedFileIDNME,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),NMERegion)
                  if NMERegion >= 0.99
                    @printf(bedFileIDENTR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),2)
                  elseif NMERegion >= 0.92
                    @printf(bedFileIDENTR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),1)
                  elseif NMERegion > 0.44
                    @printf(bedFileIDENTR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),0)
                  elseif NMERegion > 0.28
                    @printf(bedFileIDENTR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),-1)
                  elseif NMERegion >=0
                    @printf(bedFileIDENTR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),-2)
                  end #end classification NME
                end #end non-nan NME


                # Do methylation level classification

                if  (numCpGsInSubReg > 1)

                  unmethProb = cLprobsRegion[1]+cLprobsRegion[2]

                  if unmethProb > 1-thresh  # monostable unmethylated
                    if cLprobsRegion[1] > 1-thresh # highly unmethylated
                      @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),-2)
                    else # partially unmethylated
                      @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),(subEndBP-1),-1)
                    end
                  elseif unmethProb < thresh # monostable methylated
                    if cLprobsRegion[4] > 1-thresh  # highly methylated
                      @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),2)
                    else #  partially methylated
                      @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),1)
                    end
                  elseif !(sum(isnan.(cLprobsRegion))>0) # Variable
                    @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),0)  # print out 0 on METH track

                    pratio1 = cLprobsRegion[1]/unmethProb
                    pratio2 = cLprobsRegion[4]/(1-unmethProb)
                    if (pratio1 <= thresh) && (pratio2 <= thresh)
                      #mixed
                      @printf(bedFileIDVAR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),1)
                    elseif (pratio1 < 1-thresh) && (pratio2 < 1-thresh)
                      #highly mixed
                      @printf(bedFileIDVAR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),2)
                    elseif (pratio1 >= 1-thresh) && (pratio2 < 1-thresh)
                      #bistable
                      @printf(bedFileIDVAR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),3)
                    end # if VAR classification
                  end # if methylated/unmethylated/VAR classification
                else # Only one CpG site
                  unmethProb = cLprobsRegion[1]+cLprobsRegion[2]
                  if unmethProb > 1-thresh # unmethylated
                    @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),-1)
                  elseif unmethProb < thresh # methylated
                    @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),1)
                  elseif !(sum(isnan.(cLprobsRegion))>0) # mixed
                    @printf(bedFileIDMETH,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),0)
                    @printf(bedFileIDVAR,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),1)
                  end
                end # end min num of cpg for methylation classification
              end #subregion modeled
            end # loop over subregions
          end # end data present in region
        #catch
        #  println(string("Error at Chr", chr_num_str, ":bp", startBP, "-", endBP ))
        #end
      end #loop over regions
    #catch
    #  println(string("Error at Chr ", chr_num_str))
    #end
  end # loop over chromosomes

  #
  # Close all files currently open for writing
  #

  close(bedFileIDMML)
  close(bedFileIDMETH)
  close(bedFileIDVAR)
  close(bedFileIDENTR)
  close(bedFileIDNME)
  close(bedFileIDMML)
  MSIflag && close(bedFileIDMSI)
  ESIflag && close(bedFileIDESI)
  MCflag && close(bedFileIDTURN)
  MCflag && close(bedFileIDCAP)
  MCflag && close(bedFileIDRDE)
  return 0
end
