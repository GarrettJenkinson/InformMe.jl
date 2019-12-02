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
 This function makes BED files for the differential version of the
 methylation analysis results obtained by means of MethAnalysisForChr.m
 applied on two dinstict phenotypes.

 USAGE (default):

 makeBedsForDiffMethAnalysis(prefix_1,prefix_2,analysis_path_1,...
    analysis_path_2,reference_path)

 MANDATORY INPUTS:

 prefix_X
               Strings with the first string specifying
               the name of the first phenotype and the second string
               specifying the name of the second phenotype used for
               differential methylation analysis. Both phenotypes
               must have already been analyzed with methAnalysisForChr.m.

 analysis_path_X
               A string that specifies the path of the directory in which
               the methylation analysis results obtained by
               MethAnalysisForChr.m are stored.
        Default: "\$INTERMEDIATE"

 reference_path
               A string that specifies the path to the directory that
               contains the results of analysis of the reference genome
               performed by FastaToCpG.m as well as the results of
               methylation calling performed by matrixFromBam.m.
        Default: "\$REFGENEDIR"

 OPTIONAL INPUTS:

 minChrNum
               A number specifying the starting chromosome that will be
               included in the BED files.
               Default value: 1

 maxChrNum
               A number specifying the last chromosome that will be
               included in the outut BED files. Must be
               maxChrNum >= minChrNum.
               Default value: 22

 outdir
               A string that specifies the path of the directory in which
               the output BED files are written.
               Default value "\$INTERMEDIATE"

 MSIflag
               Flag that determines whether this function performs
               computation of the methylation sensitivity index (MSI).
               0: no MSI computation.
               1: allow MSI computation.
               Default value: 0

 ESIflag
               Flag that determines whether this function performs
               computation of the entropic sensitivity index (ESI).
               0: no ESI computation.
               1: allow ESI computation.
               Default value: 0

 MCflag
               Flag that determines whether this function performs
               computation of turnover ratios, CpG entropies, capacities,
               and relative dissipated energies of methylation
               channels (MCs).
               0: no MC computations.
               1: allow MC computations.
               Default value: 0

 regionSize
               The size of the genomic regions used for parameter
               estimation (in number of base pairs).
               Default value: 3000

 subregionSize
               The size of the subregions of a genomic region used
               for methylation analysis (in number of base pairs).
               The ratio regionSize/subregionSize must be an integer.
               Default value: 150

 minNumCpG     The minimum number of CpG sites within an analysis
               subregion required for performing full methylation-based
               differential classification.
               Default value: 2

 thresh

               A scalar used as a threshold in methylation-based
               differential classification.
               Default value: 0.55

 threshDMU
               A 1x6 vector containing threshold values used for
               methylation-based differential classification.
               Default value: [-1,-0.55,-0.1,0.1,0.55,1]

 threshDEU
               A 1x8 vector containing threshold values used for
               entropy-based differential classification.
               Default value: [-1,-0.5,-0.3,-0.05,0.05,0.3,0.5,1]

 The default values of regionSize, subregionSize, minNumCpG, thresh,
 threshDMU, and threshDEU should only be changed by an expert with
 a detailed understanding of the code and the methods used.

"""
function diffMethAnalysisToBed(phenoName_1,phenoName_2,analysis_path_1,analysis_path_2,reference_path,
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

  #
  # Manual checks/corrections of inputs
  #

  if analysis_path_1[end:end]!="/"
      analysis_path_1=string(analysis_path_1,"/")
  end
  if analysis_path_2[end:end]!="/"
      analysis_path_2=string(analysis_path_2,"/")
  end
  if reference_path[end:end]!="/"
      reference_path=string(reference_path, "/")
  end
  if outdir[end:end]!="/"
      outdir=string(outdir,"/")
  end

  # Hard coded options

  #MML
  hardCodedOptionsMML = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dMML-"

  #DMU
  hardCodedOptionsDMU = "track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-3.0:3.0 name=DMU-"

  #dNME
  hardCodedOptionsNME = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dNME-"

  #DEU
  hardCodedOptionsDEU = "track type=bedGraph visibility=dense windowingFunction=mean autoScale=off viewLimits=-3.0:3.0 name=DEU-"

  #JSD
  hardCodedOptionsJS = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=0.0:1.0 name=JSD-"

  #MSI
  MSIflag && (hardCodedOptionsMSI = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dMSI-")

  #ESI
  ESIflag && (hardCodedOptionsESI = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dESI-")

  #dCAP
  MCflag && (hardCodedOptionsCAP = "track type=bedGraph visibility=full windowingFunction=mean autoScale=off viewLimits=-1.0:1.0 name=dCAP-")

  #dRDE
  MCflag && (hardCodedOptionsRDE = "track type=bedGraph visibility=full windowingFunction=mean autoScale=on name=dRDE-")


  # Open .bed files to write output

  #check that bedfile folder exists, create it if not
  if !isdir(outdir)
    mkdir(outdir)
  end

 # if minChrNum > maxChrNum
 #   println("error: minChrNum must be <= maxChrNum")
 #   return 1
 # end

  # dMML
  bedFileNameMML = string(outdir, "dMML-", phenoName_1, "-VS-", phenoName_2, ".bed")
  bedFileIDMML   = open(bedFileNameMML,"w")

  # DMU
  bedFileNameDMU = string(outdir, "DMU-", phenoName_1, "-VS-", phenoName_2, ".bed")
  bedFileIDDMU   = open(bedFileNameDMU,"w")

  # dNME
  bedFileNameNME = string(outdir, "dNME-", phenoName_1, "-VS-", phenoName_2, ".bed")
  bedFileIDNME   = open(bedFileNameNME,"w")

  # DEU
  bedFileNameDEU = string(outdir, "DEU-", phenoName_1, "-VS-", phenoName_2, ".bed")
  bedFileIDDEU   = open(bedFileNameDEU,"w")

  # JSD
  bedFileNameJS = string(outdir, "JSD-", phenoName_1, "-VS-", phenoName_2, ".bed")
  bedFileIDJS   = open(bedFileNameJS,"w")

  # MSI
  MSIflag && (bedFileNameMSI = string(outdir, "dMSI-", phenoName_1, "-VS-", phenoName_2, ".bed"))
  MSIflag && (bedFileIDMSI = open(bedFileNameMSI,"w"))

  # ESI
  ESIflag && (bedFileNameESI = string(outdir, "dESI-", phenoName_1, "-VS-", phenoName_2, ".bed"))
  ESIflag && (bedFileIDESI = open(bedFileNameESI,"w"))

  # dCAP
  MCflag && (bedFileNameCAP = string(outdir, "dCAP-", phenoName_1, "-VS-", phenoName_2, ".bed"))
  MCflag && (bedFileIDCAP = open(bedFileNameCAP,"w"))

  # dRDE
  MCflag && (bedFileNameRDE = string(outdir, "dRDE-", phenoName_1, "-VS-", phenoName_2, ".bed"))
  MCflag && (bedFileIDRDE = open(bedFileNameRDE,"w"))


  # write header Info

  # dMML
  @printf(bedFileIDMML,"%s%s-VS-%s\n",hardCodedOptionsMML,phenoName_1,phenoName_2)

  # DMU
  @printf(bedFileIDDMU,"%s%s-VS-%s\n",hardCodedOptionsDMU,phenoName_1,phenoName_2)

  # dNME
  @printf(bedFileIDNME,"%s%s-VS-%s\n",hardCodedOptionsNME,phenoName_1,phenoName_2)

  # DEU
  @printf(bedFileIDDEU,"%s%s-VS-%s\n",hardCodedOptionsDEU,phenoName_1,phenoName_2)

  # JSD
  @printf(bedFileIDJS,"%s%s-VS-%s\n",hardCodedOptionsJS,phenoName_1,phenoName_2)

  # MSI
  MSIflag && @printf(bedFileIDMSI,"%s%s-VS-%s\n",hardCodedOptionsMSI,phenoName_1,phenoName_2)

  # ESI
  ESIflag && @printf(bedFileIDESI,"%s%s-VS-%s\n",hardCodedOptionsESI,phenoName_1,phenoName_2)

  # dCAP
  MCflag && @printf(bedFileIDCAP,"%s%s-VS-%s\n",hardCodedOptionsCAP,phenoName_1,phenoName_2)

  # RDE
  MCflag && @printf(bedFileIDRDE,"%s%s-VS-%s\n",hardCodedOptionsRDE,phenoName_1,phenoName_2)


  # Do the main loop over chromosomes
  for chr_str in chrs
    try
      if (length(chr_str)>3) && (chr_str(1:3)=="chr")
        chr_num_str = chr_str[4:end]
      else
        chr_num_str = chr_str
        chr_str = string("chr",chr_num_str)
      end

      # Find last CpG site on chromosome (to determine how many regions to loop through)

      # Load CpG site data on chromosome (to determine how many regions to loop through, etc.)
      CpGdata = string(reference_path, "CpGlocationChr", chr_num_str, ".jld2")
      finalCpGloc = FileIO.load(CpGdata,"finalCpGloc")
      CpGlocation = FileIO.load(CpGdata,"CpGlocation")

      #
      # find all start base-pairs for regions to be modeled on this chromsome
      #
      allStartBPs = 1:regionSize:finalCpGloc # only need to go to the last CpG site, not last BP

      #
      # Load relevant results files
      #

      analysFile1  = string(analysis_path_1, chr_str, "/", phenoName_1, "_analysis.jld2")
      mapObjAnaly1 = FileIO.load(analysFile1,"mapObjAnaly")
      analysFile2  = string(analysis_path_2, chr_str, "/", phenoName_2, "_analysis.jld2")
      mapObjAnaly2 = FileIO.load(analysFile2,"mapObjAnaly")

      ###########################################################################
      # Loop through all regions on chromosome
      ###########################################################################
      #

      for regionNum = 1:length(allStartBPs)
        try
          #
          # find path of this region
          #
          startBP = allStartBPs[regionNum]
          endBP   = startBP+regionSize-1

          locationPathName = string(chr_str, "/bp", startBP, "-", endBP)

          if haskey(mapObjAnaly1,locationPathName) && haskey(mapObjAnaly2,locationPathName)
            #
            # load structs for region
            #
            regStruct1  = mapObjAnaly1[locationPathName]
            sparseRows1 = regStruct1.sparseRows
            sparseCols1 = regStruct1.sparseCols
            sparseVals1 = regStruct1.sparseVals
            sparseNrow1 = regStruct1.sparseNrow
            sparseNcol1 = regStruct1.sparseNcol
            fullLProbs1 = SparseArrays.sparse(sparseRows1,sparseCols1,
                                              sparseVals1,sparseNrow1,
                                              sparseNcol1)
            isModeled1  = regStruct1.isModeled
            Ncg1        = regStruct1.Ncg
            NME1        = regStruct1.NME
            if ESIflag
              ESI1 = regStruct1.ESI
              if isempty(ESI1)
                  close(bedFileIDMML)
                  close(bedFileIDDMU)
                  close(bedFileIDNME)
                  close(bedFileIDDEU)
                  close(bedFileIDJS)
                  close(bedFileIDESI)
                  MSIflag && close(bedFileIDMSI)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: ESIflag = true, but the previous methylation analysis step for phenoName_1 was run with ESIflag = false.\nRerun this julia function with ESIflag = false or rerun the previous methylation analysis step with ESIflag = true after deleting $analysFile1.\n")
                return 1
              end
            end

            if MSIflag
              MSI1 = regStruct1.MSI
              if isempty(MSI1)
                  close(bedFileIDMML)
                  close(bedFileIDDMU)
                  close(bedFileIDNME)
                  close(bedFileIDDEU)
                  close(bedFileIDJS)
                  close(bedFileIDMSI)
                  ESIflag && close(bedFileIDESI)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: MSIflag = true, but the previous methylation analysis step for phenoName_1 was run with MSIflag = false.\nRerun this julia function with MSIflag = false or rerun the previous methylation analysis step with MSIflag = true after deleting $analysFile1.\n")
                return 1
              end
            end

            if MCflag
              CAP1 = regStruct1.CAP
              RDE1 = regStruct1.RDE
              if isempty(CAP1)||isempty(RDE1)
                close(bedFileIDMML)
                close(bedFileIDDMU)
                close(bedFileIDNME)
                close(bedFileIDDEU)
                close(bedFileIDJS)
                ESIflag && close(bedFileIDESI)
                MSIflag && close(bedFileIDMSI)
                close(bedFileIDCAP)
                close(bedFileIDRDE)
                println("Error: MCflag = true, but the previous methylation analysis step in phenoName_1 was run with MCflag = false.\nRerun this julia function with MCflag = false or rerun the previous methylation analysis step with MCflag = true after deleting $analysFile1.\n")
                return 1
              end
            end

            regStruct2   = mapObjAnaly2[locationPathName]
            sparseRows2 = regStruct2.sparseRows
            sparseCols2 = regStruct2.sparseCols
            sparseVals2 = regStruct2.sparseVals
            sparseNrow2 = regStruct2.sparseNrow
            sparseNcol2 = regStruct2.sparseNcol
            fullLProbs2 = SparseArrays.sparse(sparseRows2,sparseCols2,
                                              sparseVals2,sparseNrow2,
                                              sparseNcol2)
            isModeled2 = regStruct2.isModeled
            Ncg2       = regStruct2.Ncg
            NME2       = regStruct2.NME

            if MSIflag
              MSI2 = regStruct2.MSI
              if isempty(MSI2)
                  close(bedFileIDMML)
                  close(bedFileIDDMU)
                  close(bedFileIDNME)
                  close(bedFileIDDEU)
                  close(bedFileIDJS)
                  close(bedFileIDMSI)
                  ESIflag && close(bedFileIDESI)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                  println("Error: MSIflag = true, but the previous methylation analysis step for phenoName_1 was run with MSIflag = false.\nRerun this julia function with MSIflag = false or rerun the previous methylation analysis step with MSIflag = true after deleting $analysFile2.\n")
                return 1
              end
            end

            if ESIflag
              ESI2 = regStruct2.ESI
              if isempty(ESI2)
                  close(bedFileIDMML)
                  close(bedFileIDDMU)
                  close(bedFileIDNME)
                  close(bedFileIDDEU)
                  close(bedFileIDJS)
                  close(bedFileIDESI)
                  MSIflag && close(bedFileIDMSI)
                  MCflag && close(bedFileIDCAP)
                  MCflag && close(bedFileIDRDE)
                println("Error: ESIflag = true, but the previous methylation analysis step for phenoName_2 was run with ESIflag = false.\nRerun this julia function with ESIflag = false or rerun the previous methylation analysis step with ESIflag = true after deleting $analysFile2.\n")
                return 1
              end
            end

            if MCflag
              CAP2 = regStruct2.CAP
              RDE2 = regStruct2.RDE
              if isempty(CAP2)||isempty(RDE2)
                close(bedFileIDMML)
                close(bedFileIDDMU)
                close(bedFileIDNME)
                close(bedFileIDDEU)
                close(bedFileIDJS)
                MSIflag && close(bedFileIDMSI)
                ESIflag && close(bedFileIDESI)
                close(bedFileIDCAP)
                close(bedFileIDRDE)
                println("Error: MCflag = true, but the previous methylation analysis step in phenoName_2 was run with MCflag = false.\nRerun this julia function with MCflag = false or rerun the previous methylation analysis step with MCflag = true after deleting $analysFile2.\n")
                return 1
              end
            end



            if sum(Ncg1.!=Ncg2)>0
              println(string("ERROR: N1 does not equal N2 at region: ", locationPathName))
              continue
            end

            # Compute regional quantities
            lower_index,upper_index = findSortedIndices(CpGlocation,startBP,endBP)
            CpGlocInReg = CpGlocation[lower_index:upper_index]
            numCpGinReg = length(CpGlocInReg)

            subRegCount = 0
            for offset = 0:subregionSize:(regionSize-1)
              subRegCount = subRegCount + 1
              subStartBP = startBP+offset
              subEndBP   = min(subStartBP + subregionSize - 1,finalCpGloc)

              numCpGsInSubReg = Ncg1[subRegCount]

              if numCpGsInSubReg>0 && isModeled1[subRegCount] && isModeled2[subRegCount] # region modeled

                lower_index,upper_index = findSortedIndices(CpGlocInReg,subStartBP,subEndBP)
                if numCpGsInSubReg==1 && ((lower_index==1)||(upper_index==numCpGinReg)) # only 1 cpg site and its a boundary condition
                  continue #move onto next subregion
                end


                NMERegion1   = NME1[subRegCount]
                NMERegion2   = NME2[subRegCount]
                dNMERegion   = NMERegion1-NMERegion2
                if dNMERegion<Inf #no nans or Infs
                  @printf(bedFileIDNME,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),dNMERegion)
                end

                if MSIflag
                  dMSI  = MSI1[subRegCount]-MSI2[subRegCount]
                  if dMSI<Inf #no nans or infs
                    @printf(bedFileIDMSI,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),dMSI)
                  end
                end

                if ESIflag
                  dESI  = ESI1[subRegCount]-ESI2[subRegCount]
                  if dESI<Inf #no nans or infs
                    @printf(bedFileIDESI,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),dESI)
                  end
                end

                if MCflag
                  dCAP  = CAP1[subRegCount]-CAP2[subRegCount]
                  if dCAP<Inf #no nans or Infs
                    @printf(bedFileIDCAP,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),dCAP)
                  end

                  dRDE  = RDE1[subRegCount]-RDE2[subRegCount]
                  if dRDE<Inf #no nans or Infs
                    @printf(bedFileIDRDE,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),dRDE)
                  end
                end



                # Do NME classification
                if dNMERegion >= threshDEU[7] #high
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),3)
                elseif dNMERegion >= threshDEU[6] #moderate
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),2)
                elseif dNMERegion >= threshDEU[5] #weak
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),1)
                elseif dNMERegion > threshDEU[4] #iso
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),0)
                elseif dNMERegion > threshDEU[3] #weak
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),-1)
                elseif dNMERegion > threshDEU[2]   #moderate
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),-2)
                elseif dNMERegion >=threshDEU[1]   #high
                  @printf(bedFileIDDEU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),-3)
                end #end classification NME

                #
                # Compute distribution of Y1-Y2 and its mean
                #

                # TODO ########################################################
                pL1 = full(fullLProbs1[subRegCount,1:(numCpGsInSubReg+1)])
                pL2 = full(fullLProbs2[subRegCount,1:(numCpGsInSubReg+1)])
                pL1=pL1[:]#ensures column vector
                pL2=pL2[:]#ensures column vector
                # TODO ########################################################


                pD    = DSP.conv(pL1,pL2[end:-1:1])
                Dvals = -1:(1/numCpGsInSubReg):1

                #find mean
                dMML = dot(pD,Dvals)
                if dMML<Inf #no nans or Infinities
                  @printf(bedFileIDMML,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),dMML)
                end

                #
                # Compute JS distance
                #

                # JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
                # where M=(P+M)/2
                # and KL(P|M)=sum_i P_i log2(P_i/M_i)
                JSD = sqrt( sum( pL1[pL1.>0].*log2.(2 .*pL1[pL1.>0]./(pL1[pL1.>0]+pL2[pL1.>0])) )/2 +
                            sum( pL2[pL2.>0].*log2.(2 .*pL2[pL2.>0]./(pL1[pL2.>0]+pL2[pL2.>0])) )/2 )

                if JSD<Inf # ensure no nans or Infinities
                  @printf(bedFileIDJS,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                          (subEndBP-1),JSD)
                end

                #
                # Compute coarse probabilities q
                #

                qProbsRegion    = zeros(Float64,5)
                qProbsRegion[1] = sum(pD[(threshDMU[1].<=Dvals).&(Dvals.<=threshDMU[2])]) #   -1 <= Z <= -0.8
                qProbsRegion[2] = sum(pD[(threshDMU[2].<Dvals).&(Dvals.<=threshDMU[3])])  # -0.8 <  Z <= -0.1
                qProbsRegion[3] = sum(pD[(threshDMU[3].<Dvals).&(Dvals.<threshDMU[4])])   # -0.1 <  Z <   0.1
                qProbsRegion[4] = sum(pD[(threshDMU[4].<=Dvals).&(Dvals.<threshDMU[5])])  #  0.1 <= Z <   0.8
                qProbsRegion[5] = sum(pD[(threshDMU[5].<=Dvals).&(Dvals.<=threshDMU[6])]) #  0.8 <= Z <=  1

                if numCpGsInSubReg>=minNumCpG # enough to do normal classification

                  if qProbsRegion[3] > thresh
                    #ISO-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),0)

                  elseif (qProbsRegion[1]+qProbsRegion[2])/(1-qProbsRegion[3]) > thresh

                    if qProbsRegion[1] > thresh
                      #Ph-1 STRONGLY HYPO-METHYLATED
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),-3)

                    elseif (qProbsRegion[1]+qProbsRegion[2]) > thresh
                      #Ph-1 MODERATELY HYPO-METHYLATED
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),-2)

                    else
                      #Ph-1 is weakly hypomethylated
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),-1)

                    end

                  elseif (qProbsRegion[4]+qProbsRegion[5])/(1-qProbsRegion[3]) > thresh

                    if qProbsRegion[5] > thresh
                      #Ph-1 STRONGLY HYPER-METHYLATED
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),3)

                    elseif (qProbsRegion[4]+qProbsRegion[5]) > thresh
                      #Ph-1 MODERATELY HYPER-METHYLATED
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),2)

                    else
                      #Ph-1 weakly hypermethylated
                      @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                              (subEndBP-1),1)
                    end

                  end

               else # for less than minNumCpG sites. no "highly" classification allowed
                  if qProbsRegion[3] > thresh
                    #ISO-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),0)
                  elseif (qProbsRegion[4]+qProbsRegion[5]) > thresh
                    #Ph-1 HYPER-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),2)
                  elseif (qProbsRegion[1]+qProbsRegion[2]) > thresh
                    #Ph-1 HYPO-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),-2)
                  elseif (qProbsRegion[4]+qProbsRegion[5])/(1-qProbsRegion[3]) > thresh# weakly hypometh
                    #Ph-1 WEAKLY HYPER-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),1)
                  elseif (qProbsRegion[1]+qProbsRegion[2])/(1-qProbsRegion[3]) > thresh # weakly hypermeth
                    #Ph-1 WEAKLY HYPO-METHYLATED
                    @printf(bedFileIDDMU,"%s\t%u\t%u\t%f\n",chr_str,(subStartBP-1),
                            (subEndBP-1),-1)
                  end

                end # end methylation-based differential classification

              end #subregion modeled

            end # loop over subregions

          end # data in both regions
        catch
          println(string("Error at Chr", chr_num_str, ":bp", startBP, "-", endBP))
        end
      end #loop over regions
    catch
      println(string("Error at Chr ", chr_num_str))
    end
  end #loop over chromosomes

  #
  # Close all files currently open for writing
  #

  close(bedFileIDMML)
  close(bedFileIDDMU)
  close(bedFileIDNME)
  close(bedFileIDDEU)
  close(bedFileIDJS)
  MSIflag && close(bedFileIDMSI)
  ESIflag && close(bedFileIDESI)
  MCflag && close(bedFileIDCAP)
  MCflag && close(bedFileIDRDE)

  return 0
end
