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
 Compute Average Log Likelihood of Data Under Given Model

 usage (note the transpose of the dataMatrix):
 avelogLikelihood = computeAveLogLikelihood(An,Cn,dataMatrix', CpGstart, CpGend);

"""
function computeAveLogLikelihood(An::Array{Float64,1},Cn::Array{Float64,1},
                          dataMatrix::LinearAlgebra.Adjoint{Int8,Array{Int8,2}}, 
                          CpGstart::Array{<:Integer,1}, CpGend::Array{<:Integer,1})

  K = length(CpGstart)
  N = length(An)
  x_r_rPLUSs::Array{Int8,1} = zeros(Int8,N) # prealloacte for speed

  # Compute Z values
  logZ1,logZ0,logZ = computeZ(An,Cn)

  # Compute Ztilde values
  logZ1tilde,logZ0tilde,logZtilde = computeZtilde(An,Cn)

  #
  # Begin computing likelihoods for each observation
  #

  #initialize loop vars
  logAveLikeTemp::Float64 = 0 

  @inbounds for k=1:K # each observation

    r = CpGstart[k]
    s = CpGend[k]-CpGstart[k]

    # fill in observation vector
    for CpG = r:(r+s)
      x_r_rPLUSs[CpG] = dataMatrix[CpG,k] # dataMatrix is transposed before input, so stored column-wise for speedy access
    end
    
    logAveLikeTemp += calcMargProb(r,s,x_r_rPLUSs[r:(r+s)],logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn)
  end # end loop over all reads in matrix

  return (logAveLikeTemp/K)
end
