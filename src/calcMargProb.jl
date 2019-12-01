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
  Compute the Marginal Probability
 
  usage:
  logMargProb = calcMargProb(r,s,x_r_rPLUSs,logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,An,Cn)
 
"""
function calcMargProb(r::Integer,s::Integer,x_r_rPLUSs::Array{Int8,1},
              logZ1::Array{BigFloat,1},logZ0::Array{BigFloat,1},logZ::BigFloat,
              logZ1tilde::Array{BigFloat,1},logZ0tilde::Array{BigFloat,1},
              An::Array{Float64,1},Cn::Array{Float64,1}) 

  #
  # Set to log(1/Z)
  #

  logMargProb::Float64 = -logZ

  #
  # add log[\tilde{Z}_{r}(x_r)]
  #

  if x_r_rPLUSs[1]>0 # X(r)=1
    logMargProb += logZ1tilde[r]
  else # X(r)=0
    logMargProb += logZ0tilde[r]
  end


  #
  # add Z_{r+s}(x_{r+s})
  #

  if x_r_rPLUSs[s+1]>0 # X(r+s)=1
    logMargProb += logZ1[r+s]
  else                 # X(r+s)=0
    logMargProb += logZ0[r+s]
  end

  #
  # add by \sum_{n=r}^{r+s-1} log[ \phi(x_n,x_{n+1}) ]
  #

  if s>0 # otherwise this product is empty
    if r>1 # no boundary \phi term needed

      @inbounds for iter = 1:s
        n=iter+r-1
        logMargProb += ( An[n+1] + Cn[n]*(2*x_r_rPLUSs[iter]-1) )*(2*x_r_rPLUSs[iter+1]-1)
      end

    else # must compute boundary \phi term since r=1

      # n=1;
      logMargProb +=  (2*x_r_rPLUSs[1]-1)*An[1] +(2*x_r_rPLUSs[2]-1)*( An[2] +Cn[1]*(2*x_r_rPLUSs[1]-1) )

      @inbounds for iter2=2:s[1]
        logMargProb += ( An[iter2+1]+Cn[iter2]*(2*x_r_rPLUSs[iter2]-1) )*(2*x_r_rPLUSs[iter2+1]-1)
      end

    end
  end

  return logMargProb

end
