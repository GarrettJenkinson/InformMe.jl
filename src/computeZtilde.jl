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
 Backwards Recursive Computation of Partition Function

 usage:
 logZ1tilde,logZ0tilde,logZtilde = computeZtilde(An,Cn)

"""
function computeZtilde(An::Array{Float64,1},Cn::Array{Float64,1})
  #
  # Setup Multiple Precision Code
  #

  ## Required precision of computations in decimal digits
  # digitsPrec = 200
  ## convert to bits of precision
  # bitsPrec = convert(Int,ceil(digitsPrec/log10(2)))
  ## do all computations with bitsPrec precision
  #setprecision(bitsPrec) 
  setprecision(664) #hardcoded for speed  

  N = length(An)

  #
  #Initialize  output vars
  #
  logZ1tilde::Array{BigFloat,1} = zeros(BigFloat,N)
  logZ0tilde::Array{BigFloat,1} = zeros(BigFloat,N)
  #logZtilde  = zero(BigFloat)

  #
  # Begin computation of logZ's
  #

  #
  # Calculate first two boundary values
  #
# takend care of by initialization
#  logZ1tilde[1] = zero(BigFloat);
#  logZ0tilde[1] = zero(BigFloat);

  logZ1tilde[2] = log( exp(-An[1]+An[2]-Cn[1] )
                     + exp( An[1]+An[2]+Cn[1] ) )

  logZ0tilde[2] = log( exp(-An[1]-An[2]+Cn[1] )
                     + exp( An[1]-An[2]-Cn[1]) )

  #
  # Recurse through non-boundary values
  #
  @inbounds for n=3:N

      #
      # \tilde{Z}_n(x_n) = \sum_{x_{n-1}=0}^1 \phi_{n-1}(x_{n-1},x_n) \tilde{Z}_{n-1}(x_{n-1})
      #
      # This implies:
      #
      # log(\tilde{Z}_n(x_n)) = log( exp( log(\phi_{n-1}(0,x_n))+log(Z_{n-1}(0) ) ...
      #                             +exp( log(\phi_{n-1}(1,x_n))+log(Z_{n-1}(1) ))
      #

      logZ1tilde[n] = log( exp( An[n]-Cn[n-1]+logZ0tilde[n-1] )
                         + exp( An[n]+Cn[n-1]+logZ1tilde[n-1] ) )


      logZ0tilde[n] = log( exp(-An[n]+Cn[n-1]+logZ0tilde[n-1] )
                         + exp(-An[n]-Cn[n-1]+logZ1tilde[n-1] ) )

  end

  #
  # Compute Log Partition Function = log(Z_N(0)+Z_N(1))
  #

  logZtilde::BigFloat = log( exp(logZ0tilde[N]) + exp(logZ1tilde[N]) )

  return logZ1tilde,logZ0tilde,logZtilde

end
