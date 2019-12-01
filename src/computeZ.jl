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
  Recursive Computation of Partition Function

 usage:
 logZ1,logZ0,logZ=computeZ(An,Cn)

"""
function computeZ(An::Array{Float64,1},Cn::Array{Float64,1})
  #
  # Setup Multiple Precision Code
  #

  ## Required precision of computations in decimal digits
  #digitsPrec = 200
  ## convert to bits of precision
  #bitsPrec = convert(Int,ceil(digitsPrec/log10(2)))
  ## do all computations with bitsPrec precision
  #setprecision(bitsPrec) 
  setprecision(664) #hardcode for speed  

  #
  # make the inputs bigfloats for computation
  #
  N = length(An)

  #
  #Initialize  output vars
  #
  logZ1::Array{BigFloat,1} = zeros(BigFloat,N)
  logZ0::Array{BigFloat,1} = zeros(BigFloat,N)
#  logZ = zero(BigFloat)

  #
  # Begin computation of logZ's
  #

  #Calculate first boundary values 
  #initialization took care of it
  #logZ1[N] = zero(BigFloat)
  #logZ0[N] = zero(BigFloat)

  #Recurse through non-boundary values
  @inbounds for n = (N-1):-1:2 #(int n=N-1;n>1;n--){
    #log(Z_n(x_n)) = log( exp( log(\phi_n(x_n,0))+log(Z_{n+1}(0) ) ...
    #                     +exp( log(\phi_n(x_n,1))+log(Z_{n+1}(1) )) ;
    logZ1[n] = log(  exp(-An[n+1]-Cn[n]+logZ0[n+1])
                   + exp( An[n+1]+Cn[n]+logZ1[n+1])  )

    logZ0[n] = log( exp(-An[n+1]+Cn[n]+logZ0[n+1])
                   + exp( An[n+1]-Cn[n]+logZ1[n+1])  )
  end

  #Calculate last boundary values

  logZ1[1] = log( exp( An[1]-An[2]-Cn[1]+logZ0[2])
                 + exp( An[1]+An[2]+Cn[1]+logZ1[2]) )

  logZ0[1] = log( exp(-An[1]-An[2]+Cn[1]+logZ0[2])
                 + exp(-An[1]+An[2]-Cn[1]+logZ1[2]) )


  #Compute log partition function = log(Z_1(0)+Z_1(1))
  logZ::BigFloat = log( exp(logZ0[1]) + exp(logZ1[1]) )

  return logZ1,logZ0,logZ

end
