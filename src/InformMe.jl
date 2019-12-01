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
__precompile__(true)

module InformMe

#
##################################################################
# Load necessary namespaces/packages
##################################################################
#
#using ArgParse
using Distributed
import LinearAlgebra.dot
import LinearAlgebra.Adjoint
using LinearAlgebra
import Random.seed!
using QuadDIRECT
using BlackBoxOptim
using SparseArrays
using JLD2
using FileIO
using BioSequences
import Printf.@printf
import DSP.conv
import Statistics.mean

# necessary to extend isempty to MethToBits-defined-types
import Base.isempty

#
##################################################################
# Export informME functions and datatypes
##################################################################
#

export methAnalysisForChr,
#  computeAveLogLikelihood,
#  methAnalysisForRegion,
#  computeMCtransProbs,
  estParamsForChr,
#  computeLstats,
#  estParamsParallelTask,
#  computeZ,
#  estimateParams,
#  computeZtilde,
  makeBedsForMethAnalysis,
#  exactSampling,
  diffMethAnalysisToBed,
#  findSortedIndices,
#  h_func,
#  calcMargProb,
#  maxent,
#  computeAnCn,
#  processMatrix,
#  analysisParallelTask,
#  nDensity,
  fastaToCpG,
#  uniqueOct,
#  matrixFromReads,
  matrixFromBam,
#  parTaskMatrix,
#  isempty,
  EstRegStruct,
  processedBAM,
  AnalyRegStruct,
  convertBAMtoBits

#
##################################################################
# Define informME datatypes and associated methods
##################################################################
#

struct processedBAM
  # Fields
  observedMatrix::Array{Int8,2}
  #CpGlocInRegion::Array{Int64,1}

  # internal constructors
  #processedBAM(observedMatrix,CpGlocInRegion) = new(observedMatrix,CpGlocInRegion)
  # processedBAM() = new(hcat(-2, 2),[0])
  processedBAM(observedMatrix) = new(observedMatrix)
  processedBAM() = new(hcat(-2, 2))
end

# extend isempty to new type
isempty(p::processedBAM) = (p.observedMatrix[1,1] == -2)


struct EstRegStruct
  # Fields
  DistInRegion::Array{Int64,1}
  densityInRegion::Array{Float64,1}
  dataMat::Array{Int8,2}
  thetabest::Array{Float64,1}
  p0::Float64
  transProbs::Array{Float64,2}
  margProbs::Array{Float64,1}
  logZ::BigFloat

  #internal constructors
  EstRegStruct(DistInRegion,
               densityInRegion,
               dataMat,
               thetabest,
               p0,
               transProbs,
               margProbs,
               logZ) =
           new(DistInRegion,
               densityInRegion,
               dataMat,
               thetabest,
               p0,
               transProbs,
               margProbs,
               logZ)
  EstRegStruct() =
    new([-1], [0], zeros(Float64,1,1), [0], 0., zeros(Float64,1,1), [0], 0)
end

# extend isempty to new type
isempty(e::EstRegStruct) = (e.DistInRegion[1] == -1)

struct AnalyRegStruct
  # Fields
  CpGlocs_local::Array{Int64,1}
  Ncg::Array{Int16,1}
  isModeled::Array{Bool,1}
  thetabest::Array{Float64,1}
  #fullLProbs::SparseArrays.SparseMatrixCSC{Float64,Int64}
  sparseRows::Array{Int64,1}
  sparseCols::Array{Int64,1}
  sparseVals::Array{Float64,1}
  sparseNrow::Int64
  sparseNcol::Int64
  cLprobs::Array{Float64,2}
  MML::Array{Float64,1}
  NME::Array{Float64,1}
  MSI::Array{Float64,1}
  ESI::Array{Float64,1}
  TURN::Array{Float64,1}
  CAP::Array{Float64,1}
  RDE::Array{Float64,1}

  # internal constructors
  AnalyRegStruct(CpGlocs_local,
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
                RDE)=
           new(CpGlocs_local,
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

  AnalyRegStruct()=
    new(Int64[],Int16[],Bool[],Float64[],
        #SparseArrays.sparse(Int64[],Int64[],Float64[],1,1),
        Int64[],Int64[],Float64[],1,1,
        zeros(Float64,1,1),Float64[],Float64[],
        Float64[],Float64[],Float64[],Float64[])
end

# extend isempty to new type
isempty(a::AnalyRegStruct) = isempty(a.Ncg)

#
##################################################################
# Load source code with all the functions
##################################################################
#

include("analysisParallelTask.jl")
include("convertBAMtoBits.jl")
include("estParamsForChr.jl")
include("estParamsParallelTask.jl")
include("estimateParams.jl") # fully tested
include("fastaToCpG.jl") # fully tested
include("matrixFromBam.jl")
include("matrixFromReads.jl") # fully tested
include("calcMargProb.jl") # fully tested
include("computeAnCn.jl") # fully tested
include("computeAveLogLikelihood.jl") # fully tested
include("computeLstats.jl") # fully tested
include("computeMCtransProbs.jl") # fully tested
include("computeZ.jl") # fully tested
include("computeZtilde.jl") # fully tested
include("diffMethAnalysisToBed.jl")
include("exactSampling.jl") # fully tested
include("findSortedIndices.jl") # fully tested
include("h_func.jl") # fully tested
include("makeBedsForMethAnalysis.jl")
include("maxent.jl") # fully tested
include("methAnalysisForRegion.jl")
include("methAnalysisForChr.jl")
include("nDensity.jl") # fully tested
include("parTaskMatrix.jl") # fully tested
include("processMatrix.jl") # fully tested
include("uniqueOct.jl") # fully tested


end # module
