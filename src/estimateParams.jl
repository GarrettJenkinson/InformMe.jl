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
 This function takes a list of BAM files (which correspond to the same
 phenotype) and estimates the parameters of the 1D Ising model that
 best fits the methylation data associated with a specific genomic
 region.

 USAGE:

 regionStruct = estimateParams(locationPathName,phenoName,...
                                 DistInRegion,densityInRegion,dataMat)

 INPUTS:

 locationPathName
               A string that specifies the relative path to the genomic
               region being modeled. For example: 'chr1/bp1-3000'.

 phenoName
               A string that specifies the name of the modeled phenotype.

 DistInRegion
               An Nx1 vector of distances to the next CpG site for each
               of the N CpG sites in the genomic region being modeled.

 densityInRegion
               An Nx1 vector of densities for each of the N CpG sites in
               the genomic region being modeled.

 dataMat
               A dxN matrix of methylation data within the genomic region
          calcMargProb(r,s,x_r_rPLUSs,logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn)     being modeled. Each row of this matrix is a single
               observation with elements taking values -1,0,1
               corresponding to the methylation status of a CpG site
               not being observed (-1), being unmethylated (0), or being
               methylated (1).

 OPTIONAL INPUTS:

 boundaryConditions
              Flag to decide if boundary conditions should be estimated
              freely in MLE.
              Default value: false

useQuad
              Flag to decide if QUADdirect should be used for optimization
              (if true) or if BlackBoxOptim should be used (if false).
              Default value: true
"""
function  estimateParams(locationPathName::String,phenoName::String,
                          DistInRegion::Array{<:Integer,1},densityInRegion::Array{Float64,1},
                          dataMat::Array{<:Integer,2};boundaryConditions=false,useQuad=true)

  #Dist::Array{Float64,1}  = DistInRegion[1:end-1] # converts to Float
  #density = deepcopy(densityInRegion)
  N       = length(densityInRegion)

  emptyStruct = EstRegStruct()
  #
  ###########################################################################
  # Check for sufficient data coverage and writeout coverage to data location
  ###########################################################################
  #

  percentCovered = sum(sum(dataMat.>-1,dims=1).>0)/N
  depthOfCov     = sum(dataMat.>-1)/N
  if (depthOfCov<2.5) || (percentCovered<(2/3))
    return  emptyStruct # insufficient data; do not build statistical model
  end

  #
  ###########################################################################
  # Process data matrix to have only contiguous reads
  ###########################################################################
  #
#  try
    newMatrix, CpGstart, CpGend = processMatrix( dataMat )
#  catch
#    println(string("Error processing data matrix from: ", locationPathName))
#    return emptyStruct
#  end

  ###########################################################################
  # Compute maximum likelihood estimator
  ###########################################################################


  # first define the objective function to be minimized
  function objFnToMinimize(theta::AbstractArray)

    # Compute model from parameters
    An,Cn = computeAnCn(densityInRegion,DistInRegion[1:(end-1)],theta)

    # Compute average log likelihood (function to maximize)
    aveLogLikelihood::Float64 = computeAveLogLikelihood(An,Cn,newMatrix', CpGstart, CpGend)
    # println("$theta,$aveLogLikelihood")
    # Objective function to minimize is the negative of the function to maximize
    return -1*aveLogLikelihood
  end

#  try

  if !useQuad
    # assumes "using BlackBoxOptim"
    
    if boundaryConditions
      optRes = bboptimize(objFnToMinimize;
                          SearchRange = [(-10.0,10.0),
                                         (-100.0,100.0),
                                         (-20.0,20.0),
                                         (-3.0,3.0),
                                         (-3.0,3.0)], 
                          Method = :dxnes, 
                          MaxTime = 10.0)
    else
      optRes = bboptimize(objFnToMinimize;
                          SearchRange = [(-10.0,10.0),
                                         (-100.0,100.0),
                                         (-20.0,20.0)], 
                          Method = :dxnes,
                          MaxTime = 10.0)
    
    end
    thetabest = best_candidate(optRes)
  else
    # assumes "using QuadDIRECT"
    if boundaryConditions
      lower = [-10.0,
              -100.0,
               -20.0,
                -3.0,
                -3.0]
      upper = [10.0,
              100.0,
               20.0,
                3.0,
                3.0]
      splits = ([-3.0,0.0,3.0],
                [-20.0,0.0,20.0],
                [-5.0,0.0,5.0],
                [-2,0.0,2],
                [-2,0.0,2])
    else
      lower = [-10.0,
              -100.0,
               -20.0]
      upper = [10.0,
              100.0,
               20.0]
      splits = ([-3.0,0.0,3.0],
                [-20.0,0.0,20.0],
                [-5.0,0.0,5.0])
    end
    nParams = length(lower)
    nQnewton = (nParams+1)*(nParams+2)*3/2
    root, x0 = analyze(objFnToMinimize, splits, lower, upper)#;
                      # maxevals=12000, print_interval=500,nquasinewton=nQnewton,rtol=0.0)
    box = minimum(root)
    fbest = value(box)
    thetabest =  position(box, x0)
  end # end useQuad
#  catch
#    println(string("Error in QuadDirect for phenotype ", phenoName, " at location: ", locationPathName))
#    return emptyStruct
#  end


#  try
    #
    ###########################################################################
    # Solve Model
    ###########################################################################
    #

    An,Cn = computeAnCn(densityInRegion,DistInRegion[1:end-1],thetabest)
    logZ1,logZ0,logZ = computeZ(An,Cn)
    logZ1tilde,logZ0tilde,logZtilde = computeZtilde(An,Cn)
    p0,transProbs = computeMCtransProbs(An,Cn,logZ1,logZ0,logZ)

    #
    # Compute 1D marginals
    #
    margProbs = zeros(Float64,N) #P(X_n=1)
    margProbs[1] = 1-p0

    for r=2:N
      #s = 0; #x_r_rPLUSs = 1;
      logMargProb = calcMargProb(r,0,Int8[1],logZ1,logZ0,logZ,
                                 logZ1tilde,logZ0tilde,An,Cn)
      margProbs[r] = exp(logMargProb)
    end

    #
    # correct numerical errors
    #
    function correctProbs(prob::Float64)
      prob > (1-eps()) && return 1-eps()
      prob < eps() && return eps()
      return prob
    end
    margProbs = correctProbs.(margProbs)
    transProbs = correctProbs.(transProbs)

#  catch
#    println(string("Error in Model Solving for phenotype ", phenoName, " at location: ", locationPathName));
#    return emptyStruct
#  end

  #
  # Store results in a structure to return
  #
#  try

    regionStruct = EstRegStruct(DistInRegion,
                                densityInRegion,
                                dataMat,
                                thetabest,
                                p0,
                                transProbs,
                                margProbs,
                                logZ)

#  catch
#    println(string("Error in creating structure for phenotype ", phenoName, " at location: ", locationPathName))
#    return emptyStruct
#  end
  return regionStruct
end
