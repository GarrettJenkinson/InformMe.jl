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
"""
 This function performs methylation analysis of a genomic region used
 for estimating the parameters of the Ising model by computing a number
 of statistical summaries of the methylation state within the region,
 including probability distributions of methylation levels, mean
 meathylation levels, and normalized methylation entropies. If desired,
 this function also computes entropic sensitivity indices, as well
 information-theoretic quantities associated with methylation channels,
 such as trunover ratios, channel capacities, and relative dissipated
 energies.

 USAGE:

 regionStruct = MethAnalysisForRegion(localEstStruct,CpGlocs_local,...
                              startBP,endBP,subregionSize,MSIflag,ESIflag,MCflag)

 INPUTS:

 localEstStruct
               A structure generated during parameter estimation within
               the genomic region.

 CpGlocs_local
               A vector of CpG locations within the genomic region.

 startBP
               The starting base pair index (1-based) along the genome for
               the genomic region.

 endBP
               The ending base pair index (1-based) along the genome for
               the genomic region.

 subregionSize
               A scalar that specifies the number of base pairs within
               subregions of the genomic region determining the
               resolution of methylation analysis. The analysis subregions
               must be of the same length (in base pairs) and chosen
               so that the genomic region is partitioned into an integer
               number of nonoverlapping subregions.

 MSIflag
               Flag that determines whether this function performs
               computation of the methylation sensitivity index (MSI).
               false: no MSI computation.
               true: allow MSI computation.
               Default value: false

 ESIflag
               Flag that determines whether this function performs
               computation of the entropic sensitivity index (ESI).
               false: no ESI computation.
               true: allow ESI computation.
               Default value: false

 MCflag
               Flag that determines whether this function performs
               computation of turnover ratios, CpG entropies, capacities,
               and relative dissipated energies of methylation
               channels (MCs).
               false: no MC computations.
               true: allow MC computations.
               Default value: false

 boundaryConditions
              Flag to decide if boundary conditions should be estimated
              freely in MLE.
              Default value: false
 OUTPUT:

 regionStruct
               A structure summarizing the methylation analysis results
               containing the following information:
               o The locations of the CpG sites within the genomic region.
               o Numbers of CpG sites within the analysis subregions.
               o Which analysis subregions are modeled and which are not.
               o Estimated parameters of Ising model in genomic region
               o Methylation level probabilities in modeled subregions.
               o Coarse methylation level probabilities.
               o Mean methylation levels.
               o Normalized methylation entropies.
               o Methylation sensitivity indices (if MSIflag).
               o Entropic sensitivity indices (if ESIflag).
               o Turnover ratios (if MCflag).
               o Channel capacities (if MCflag).
               o Relative dissipated energies (if MCflag).

"""
function methAnalysisForRegion(localEstStruct,CpGlocs_local,startBP,endBP,
                        subRegionSize;MSIflag=false,ESIflag=false,MCflag=false,boundaryConditions=false)

  subRegStartBPlist = startBP:subRegionSize:endBP

  numSubRegions = length(subRegStartBPlist)

  isModeled   = zeros(Bool,numSubRegions)

  thresh      = [0,0.25,0.5,0.75,1]

  cLprobs     = zeros(Float64,numSubRegions,length(thresh)-1)
  Ncg         = zeros(Int16,numSubRegions)
  MML         = zeros(Float64,numSubRegions)
  NME         = zeros(Float64,numSubRegions)

  if MSIflag
    MSI = zeros(Float64,numSubRegions)
  else
    MSI = Float64[]
  end

  if ESIflag
    ESI = zeros(Float64,numSubRegions)
  else
    ESI = Float64[]
  end

  if MCflag
    TURN = zeros(Float64,numSubRegions)
    CAP  = zeros(Float64,numSubRegions)
    RDE  = zeros(Float64,numSubRegions)
  else
    TURN = Float64[]
    CAP  = Float64[]
    RDE  = Float64[]
  end

  # sparse matrix parameters
  sparseRows::Array{Int64,1}   = fill(0,convert(Int64,floor(subRegionSize/2)*numSubRegions))
  sparseCols::Array{Int64,1}   = fill(0,convert(Int64,floor(subRegionSize/2)*numSubRegions))
  sparseVals::Array{Float64,1} = fill(0.0,convert(Int64,floor(subRegionSize/2)*numSubRegions))
  sparseIndex = 1

  #
  ###########################################################################
  # initialize model from data structure
  ###########################################################################
  #
  thetabest       = localEstStruct.thetabest
  margProbs       = localEstStruct.margProbs
  DistInRegion    = localEstStruct.DistInRegion
  densityInRegion = localEstStruct.densityInRegion
  transProbs      = localEstStruct.transProbs

  #correct numerical errors:
  margProbs[margProbs.<0].=0.0
  margProbs[margProbs.>1].=1.0

  if ESIflag || MSIflag
    epsilon = 0.01
    N = length(CpGlocs_local)

    Ana,Cna = computeAnCnm(densityInRegion,DistInRegion[1:(end-1)],thetabest.*(1+epsilon))

    logZ1a,logZ0a,logZa = computeZ(Ana,Cna)
    logZ1tildea,logZ0tildea,logZtildea = computeZtilde(Ana,Cna)
    p0a,transProbsa = computeMCtransProbs(Ana,Cna,logZ1a,logZ0a,logZa)

    # Compute 1D marginals
    margProbsa = zeros(Float64,N) #P(X_n=1)
    margProbsa[1] = 1-p0a

    for r=2:N
      #s = 0; #x_r_rPLUSs = 1;
      logMargProba = calcMargProb(r,0,Int8[1],logZ1a,logZ0a,logZa,
                                  logZ1tildea,logZ0tildea,Ana,Cna)
      margProbsa[r] = convert(Float64,exp(logMargProba))
    end

    # correct numerical errors
    margProbsa[margProbsa.>(1-eps())]   .= 1-eps()
    margProbsa[margProbsa.<eps()]       .= eps()
    transProbsa[transProbsa.>(1-eps())] .= 1-eps()
    transProbsa[transProbsa.<eps()]     .= eps()
  end

  ###########################################################################
  # loop through regions
  ###########################################################################

  subRegCount = 0

  for subRegStartBP = subRegStartBPlist

    subRegCount = subRegCount+1
    subRegEndBP = subRegStartBP+subRegionSize-1

    # find CpG sites
    lower_index,upper_index = findSortedIndices(CpGlocs_local,subRegStartBP,subRegEndBP)
    CpGlocs                 = CpGlocs_local[lower_index:upper_index]
    Ncg[subRegCount]        = length(CpGlocs)

    # calc methylation level probabilities in subregion
    if Ncg[subRegCount]>0

      isModeled[subRegCount] = true # set flag, this region is modeled

      # Compute mean methylation level
      MML[subRegCount]        = mean(margProbs[lower_index:upper_index])

      # Compute NME
      if Ncg[subRegCount] >= 2 # 2 or more CpGs

        LProbs,LVals,NMETemp = computeLstats(margProbs[lower_index],transProbs[lower_index:(upper_index-1),:])

        if sum(isnan.(LProbs))>0
          isModeled[subRegCount] = false
          continue
        end

        NME[subRegCount]=NMETemp

        # Add probability to sparse matrix:
        # fullLProbs[subRegCount,1:length(LProbs)]=LProbs

        for indexNum=1:length(LProbs)
          sparseVals[sparseIndex]=LProbs[indexNum]
          sparseRows[sparseIndex]=subRegCount
          sparseCols[sparseIndex]=indexNum
          sparseIndex += 1
        end


        # Compute coarse L probabilities
        if !isempty(LProbs[LVals.==thresh[3]])
          correction = LProbs[LVals.==thresh[3]]/2
        else
          correction = 0.0
        end

        cLprobs[subRegCount,1] = sum(LProbs[(thresh[1].<=LVals).&(LVals.<=thresh[2])]) # 0<= L <= 0.25
        cLprobs[subRegCount,2] = sum(LProbs[(thresh[2].<LVals).&(LVals.<thresh[3])]) + correction[1] # 0.25< L <= 0.5
        cLprobs[subRegCount,3] = sum(LProbs[(thresh[3].<LVals).&(LVals.<thresh[4])]) + correction[1]  # 0.5<= L < 0.75
        cLprobs[subRegCount,4] = sum(LProbs[(thresh[4].<=LVals).&(LVals.<=thresh[5])]) # 0.75=< L <= 1

        if ESIflag || MSIflag
          LProbsa,dummy,NMEa = computeLstats(margProbsa[lower_index],transProbsa[lower_index:(upper_index-1),:])
        end

        if MSIflag
          #do sensitivty calculation
          pL1=LProbs
          pL2=LProbsa
          # JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
          # where M=(P+M)/2
          # and KL(P|M)=sum_i P_i log2(P_i/M_i)
          JSD = sqrt( sum( pL1[pL1.>0].*log2.(2 .*pL1[pL1.>0]./(pL1[pL1.>0]+pL2[pL1.>0])) )/2 +
                      sum( pL2[pL2.>0].*log2.(2 .*pL2[pL2.>0]./(pL1[pL2.>0]+pL2[pL2.>0])) )/2 )
          MSI[subRegCount] = JSD/epsilon
        end

        if ESIflag
          #do sensitivity calculation
          if sum(isnan.(NMEa))==0
            Da=abs(NMEa-NME[subRegCount])*log2(Ncg[subRegCount]+1)
            ESI[subRegCount] = Da/epsilon
          end
        end #ESIflag

      else# Ncg[subRegCount]==1

        # just use exact marginal probabilities

        cLprobs[subRegCount,1]=1-margProbs[lower_index]
        cLprobs[subRegCount,4]=margProbs[lower_index]

        NME[subRegCount] = -margProbs[lower_index]*log2(margProbs[lower_index]) -
          (1-margProbs[lower_index])*log2(1-margProbs[lower_index])

        # Add probs to sparse matrix:
        # fullLProbs(subRegCount,1:2)=[1-margProbs(lower_index),margProbs(lower_index)]

        sparseVals[sparseIndex]   = 1-margProbs[lower_index]
        sparseVals[sparseIndex+1] = margProbs[lower_index]
        sparseRows[sparseIndex]   = subRegCount
        sparseRows[sparseIndex+1] = subRegCount
        sparseCols[sparseIndex]   = 1
        sparseCols[sparseIndex+1] = 2
        sparseIndex=sparseIndex+2

        if MSIflag || ESIflag
          LProbsa = [1-margProbsa[lower_index],margProbsa[lower_index]]
          LProbs  = [1-margProbs[lower_index], margProbs[lower_index]]
        end

        if MSIflag
          #do sensitivty calculation
          pL1=LProbs
          pL2=LProbsa
          # JSDIS := sqrt(KL(P|M)/2 + KL(Q|M)/2)
          # where M=(P+M)/2
          # and KL(P|M)=sum_i P_i log2(P_i/M_i)
          JSD = sqrt( sum( pL1[pL1.>0].*log2.(2 .*pL1[pL1.>0]./(pL1[pL1.>0]+pL2[pL1.>0])) )/2 +
                      sum( pL2[pL2.>0].*log2.(2 .*pL2[pL2.>0]./(pL1[pL2.>0]+pL2[pL2.>0])) )/2 )
          MSI[subRegCount] = JSD/epsilon
        end

        if ESIflag
          # Do sensitivity calc
          NMEa = LinearAlgebra.dot(-LProbsa,log2.(LProbsa))
          Da = abs(NMEa-NME[subRegCount])
          ESI[subRegCount] = Da/epsilon
        end #ESIflag

      end # end  if Ncg(subRegCount) >=2 else Ncg(subRegCount)==1

      if MCflag
        # Compute Capacities and TURN variables
        lambdaVals  = margProbs[lower_index:upper_index]./
                      (1 .- margProbs[lower_index:upper_index]) # P(1)/P(0)

        TURN[subRegCount] = mean(log2.(lambdaVals))

        # Compute Capacity
        CAPvals = fill(0.0,size(lambdaVals))

        CAPvals[lambdaVals.>=1] .= 1 .- 0.52 .* h_func( lambdaVals[lambdaVals.>=1]./
                                                     (1 .+ lambdaVals[lambdaVals.>=1]) )./
                                             (1 .+ lambdaVals[lambdaVals.>=1])
        CAPvals[lambdaVals.<1] .= 1 .- 0.52 .* h_func( lambdaVals[lambdaVals.<1]./
                                                    (1 .+ lambdaVals[lambdaVals.<1]) ).*
                                             lambdaVals[lambdaVals.<1]./
                                             (1 .+ lambdaVals[lambdaVals.<1])
        CAP[subRegCount] = mean(CAPvals)

        # Compute average RDE in region

        RDEvals = zeros(lambdaVals)
        RDEvals[lambdaVals.<1]  .= log2.( (1 .+ lambdaVals[lambdaVals.<1])./
                                        (2 .*lambdaVals[lambdaVals.<1])    ) .+
                                  4.76
        RDEvals[lambdaVals.>=1] .= log2.( (1 .+ lambdaVals[lambdaVals.>=1])./2 ) .+
                                  4.76
        RDE[subRegCount] = mean(RDEvals)

      end #end MCflag

    end# end Ncg(subRegCount)>0 condition

  end # end loop over subregions

  # Create sparse matrix fullLProbs with rows containing the probability
  # distribution of L

  sparseIndex=sparseIndex-1
  # fullLProbs::SparseArrays.SparseMatrixCSC{Float64,Int64} =
  #           SparseArrays.sparse(sparseRows[1:sparseIndex],sparseCols[1:sparseIndex],
  #                                  sparseVals[1:sparseIndex],numSubRegions,
  #                                  max(sparseCols[1:sparseIndex]...))
  sparseRows = sparseRows[1:sparseIndex]
  sparseCols = sparseCols[1:sparseIndex]
  sparseVals = sparseVals[1:sparseIndex]
  sparseNrow = numSubRegions
  sparseNcol = max(sparseCols...)

  # Save all results into a structure

  regionStruct = AnalyRegStruct(CpGlocs_local,
                                Ncg,
                                isModeled,
                                thetabest,
                                #fullLProbs,
                                sparseRows,
                                sparseCols,
                                sparseVals,
                                sparseNrow,
                                sparseNcol,
                                cLprobs,
                                MML,
                                NME,
                                MSI,
                                ESI,
                                TURN,
                                CAP,
                                RDE)

  return regionStruct
end
