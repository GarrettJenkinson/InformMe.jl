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
 Recursive Computation of Markov Chain Transition Probabilities

 usage:
 p1,transProbs = computeMCtransProbs(An,Cn,logZ1,logZ0,logZ)

 INPUTS:
      Omega_alpha is a Nx1 vector of Omega*alpha_n parameters for n=1,...,N
      z1 is a Nx1 vector of values Z_n(x_n) for x_n=1 and n=1,...,N
      z0 is a Nx1 vector of values Z_n(x_n) for x_n=0 and n=1,...,N
      Z is the partition function

  OUTPUTS:
      p1 is the marginal probability  P(X_1=0)
      transProbs is a (N-1)x2 matrix
        -First  column has probability P(X_{n+1}=0|X_n=0) for n=1,2,...,N-1
        -Second column has probability P(X_{n+1}=0|X_n=1) for n=1,2,...,N-1

      Note the joint dist can be calculated from these probabilities:
      P(X=x) = P(X_1) prod_{n=1}^{N-1} P(X_{n+1}|X_n)

      More importantly these probabilities can be used to iteratively draw an
      exact sample from the joint distribution without MCMC.
"""
function computeMCtransProbs(An::Array{Float64,1},Cn::Array{Float64,1},
                             logZ1::Array{BigFloat,1},logZ0::Array{BigFloat,1},logZ::BigFloat)
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
  tempVar = zero(BigFloat)

  #
  #Initialize  output vars
  #
  transProbs::Array{Float64,2} = zeros(Float64,N-1,2)

  #
  # Compute boundary condition
  #

  # P(X_1=0)
  p1::Float64  = exp( logZ0[1]-logZ )

  #
  # P(x_2=0|x_1=0) = \phi_1(x_1=0,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=0)
  #

  transProbs[1,1] = exp( -An[1]-An[2]+Cn[1]+logZ0[2]-logZ0[1] )

  #
  # P(x_2=0|x_1=1) = \phi_1(x_1=1,x_{2}=0)*z_{2}(x_{2}=0)/z_1(x_1=1)
  #
  transProbs[1,2] = exp( An[1]-An[2]-Cn[1]+logZ0[2]-logZ1[1] )

  #
  # Use backwards recursion to compute transition probabilities
  #
  @inbounds for n=2:(N-1)
      #
      # P(x_{n+1}=0|x_n=0) = \phi_n(x_n=0,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=0)
      #
      transProbs[n,1] = exp( -An[n+1]+Cn[n]+logZ0[n+1]-logZ0[n] )

      #
      # P(x_{n+1}=0|x_{n}=1) = \phi_n(x_n=1,x_{n+1}=0)*z_{n+1}(x_{n+1}=0)/z_n(x_n=1)
      #
      transProbs[n,2] = exp( -An[n+1]-Cn[n]+logZ0[n+1]-logZ1[n] )
  end

  return p1,transProbs

end

