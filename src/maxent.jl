# Permission was provided by the author, Ali Mahammad-Djafari, 
# to modify and distribute this code with the informME method's
# code. The original MATLAB code has been ported by Garrett Jenkinson
# to the julia language. Contact the original author 
# directly for use outside of informME.
#
# Author's website:
# http://djafari.free.fr/index.htm
# 
# Author's paper with source code:
#Mohammad-Djafari A. (1992) A Matlab Program to Calculate the Maximum Entropy Distributions. In: Smith C.R., Erickson G.J., Neudorfer P.O. (eds) Maximum Entropy and Bayesian Methods. Fundamental Theories of Physics (An International Book Series on The Fundamental Theories of Physics: Their Clarification, Development and Application), vol 50. Springer, Dordrecht

"""
   usage:
   lambda,p,entr=maxent(mu1,x)
  
   Modified by Garrett Jenkinson to work in julia
   with additional speedup modifications
   Original script came with the following documentation:
  
  ME_DENS2
   [LAMBDA,P,ENTR]=ME_DENS2(MU,X,LAMBDA0)
   This program calculates the Lagrange Multipliers of the ME
   probability density functions p(x) from the knowledge of the
   N moment contstraints in the form:
   E{x^n}=mu(n) n=0:N with mu(0)=1.
  
   MU is a table containing the constraints MU(n),n=1:N.
   X is a table defining the range of the variation of x.
   LAMBDA0 is a table containing the first estimate of the LAMBDAs.
   (This argument is optional.)
   LAMBDA is a table containing the resulting Lagrange parameters.
   P is a table containing the resulting pdf p(x).
   ENTR is a table containing the entropy values at each
   iteration.
  
   Author: A. Mohammad-Djafari
   Date : 10-01-1991
"""
function maxent(mu1::Array{Float64,1},xrow::AbstractArray)
  maxIters = 5000
  entr     = zeros(Float64,maxIters)

  mu = cat(1.0,mu1[:],dims=1)   # copy and add mu(0)=1
  x=collect(xrow) # copy and ensure column vector
  lx=length(x) # x axis
  xmin=x[1]
  xmax=x[lx]
  dx=x[2]-x[1]

  # initialize LAMBDA
  lambda    = fill(0.0,size(mu))         # This produces a uniform
  lambda[1] = log(xmax-xmin) # distribution.

  N=length(lambda)
  #
  M=2*N-1 # Calcul de fin(x)=x.^n
  fin=zeros(Float64,length(x),M) #
  fin[:,1] .= 1.0 # fi0(x)=1
  @inbounds for n=2:M
    fin[:,n]=x.*fin[:,n-1]
  end
  #
  iter=0
  @inbounds while true # start iterations
    iter += 1
    #
    p = exp.(-(fin[:,1:N]*lambda)) # Calculate p(x)
    #
    G=zeros(Float64,M) # Calculate Gn
    for n=1:M
      G[n]=dx*sum(fin[:,n].*p)
    end
    #
    entr[iter]=dot(lambda,G[1:N]) # Calculate the entropy value
    #
    gnk=zeros(Float64,N,N) # Calculate gnk
    for i=1:N # Matrix G is a Hankel matrix
      gnk[:,i]=(-1).*G[i:N+i-1]
    end
    #
    v=mu-G[1:N] # Calculate v
    delta=gnk\v # Calculate delta
    lambda=lambda+delta # Calculate lambda
    eps=1e-6 # Stopping rules
    if sum(abs.(delta./lambda))<eps
      break
    end
    if iter>2 && abs((entr[iter]-entr[iter-1])/entr[iter])<eps
      break
    end
    if iter>maxIters
      println("Warning: Maximum number of iterations reached")
      break
    end
  end
  #
  p=exp.(-(fin[:,1:N]*lambda)) # Calculate the final p(x)
  entr=entr[1:iter]
  p = p./sum(p)
  return lambda,p,entr
end
