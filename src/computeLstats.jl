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
#
"""
 This function computes the probability distribution of the methylation 
 level L within a genomic region comprised N CpG sites as well as the 
 corresponding normalized methylation entropy.

 USAGE:

 [LProbs,LVals,NME] = computeLstats(p1,transProbs)

 INPUTS:

 p1        
           The probability Pr[X_1=1] of the first CpG site in the genomic 
           region to be methylated. 

 transProbs 
           (N-1)x2 matrix whose elements are given by:  
               transProbs(n,1) = Pr[X_{n+1}=0 | X_n=0]
               transProbs(n,2) = Pr[X_{n+1}=0 | X_n=1]
            where X_n is the methylation status of the n-th CpG 
            site within the genomic region containing N CpG sites.

 OUTPUTS:

 LProbs    
           (N+1)x1 vector containing the probabilities of the 
            methylation level L to take values 0,1/N,...,1.

 LVals     
           (N+1)x1 vector containing the methylation level 
            values 0,1/N,...,1.

 NME       The normalized methylation entropy.

"""
function computeLstats(p1::Float64,transProbs::Array{Float64,2})


  #
  # Initialize
  #

  seed!(123)      # use common random numbers (same seed will be used everytime)
  # useful in ESI for reducing variance; see Spall
  # "Introduction to Stochastic Search and Optimization" for
  # more details on CRNs.

  maxMom              = 4   # number of moments used in MaxEnt Computations
  totMCsamps          = 2^17# number of Monte Carlo samples used in estimation
  threshForAnalytical = 18 #log2.(totMCsamps).+1
  N                   = size(transProbs,1)+1 # number of CpG sites

  transProbs[transProbs.<0] .= 0 # fix numerical errors
  transProbs[transProbs.>1] .= 1 # fix numerical errors

  #
  # Compute Yvals
  #
  yVals = 0.0:(1/N):1.0  # linspace(0,1,N+1)

  #
  # Compute Yprobs
  #
  if N<=threshForAnalytical # compute analytically

    yProbs = zeros(Float64,N+1)

    @inbounds for xPatNum=0:((2^N)-1)
      #
      # compute probability of this x pattern
      #
      xBits = string(xPatNum,base=2,pad=N) #bin(xPatNum,N) # convert to string with N bit representation of xPatNum

      sumOfX = 0

      if xBits[1]=='1' # if first bit is a 1
        xProb   = p1
        prevBit = '1'
        sumOfX += 1
      elseif xBits[1]=='0'
        xProb   = 1-p1
        prevBit = '0'
      else
        println("Error in computing Y probs")
      end

      for n=2:N
        #
        # Update probability to transition to next bit
        #
        currBit = xBits[n]

        if (currBit=='0')&&(prevBit=='0')
          xProb  = xProb*transProbs[n-1,1]    # Pr[x_{n}=0|x_{n-1}=0] Pr[x_{n-1}=0,..x_{1}]
        elseif (currBit=='0')&&(prevBit=='1')
          xProb  = xProb*transProbs[n-1,2]    # Pr[x_{n}=0|x_{n-1}=1] Pr[x_{n-1}=1,..x_{1}]
        elseif (currBit=='1')&&(prevBit=='0')
          xProb  = xProb*(1-transProbs[n-1,1])# Pr[x_{n}=1|x_{n-1}=0] Pr[x_{n-1}=0,..x_{1}]
          sumOfX = sumOfX+1
        elseif (currBit=='1')&&(prevBit=='1')
          xProb  = xProb*(1-transProbs[n-1,2])# Pr[x_{n}=1|x_{n-1}=1] Pr[x_{n-1}=1,..x_{1}]
          sumOfX = sumOfX+1
        else
          println("Error in computing Y probs")
        end

        # update previous bit for next step in loop
        prevBit=currBit
      end # end loop over N bits

      yProbs[sumOfX+1] = yProbs[sumOfX+1] + xProb

    end # end loop over patterns of X

    # correct for numerical problems:
    yProbs[yProbs.<0] .= 0.0
    yProbs = yProbs./sum(yProbs)

  else # N>18, so compute via MaxEnt and Monte Carlo
    #
    # Calculate moments using M monte carlo samples
    #
    moments = zeros(Float64,maxMom)
    @inbounds for m=1:totMCsamps

      Xsamp    = exactSampling(p1,transProbs)
      fracMeth = mean(Xsamp)

      for momNum=1:maxMom
        moments[momNum]=moments[momNum]+(fracMeth^momNum)
      end

    end# end monte carlo sample loop

    moments = moments./totMCsamps

    #
    # Compute the maxEnt model
    #

    ignore1,p,ignore2 = maxent(moments,yVals)
    p[p.<0] .= 0
    yProbs  = p./sum(p)

  end # end decision to compute analytically or via MaxEnt

  #
  # NME of Y
  #
  yProbsNonZero = yProbs[yProbs.>0]
  entropy       = dot(yProbsNonZero,(-1).*log2.(yProbsNonZero))

  NMEy = entropy/log2(N+1)
  if NMEy<0
    NMEy=0
  elseif NMEy>1
    NMEy=1
  end

  return yProbs,yVals,NMEy
end # end of function
