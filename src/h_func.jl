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
 This function computes the (log2-based) entropies of a collection 
 X(q,r), q = 1,2,...,Q, r = 1,2,...,R of QxR binary random variables 
 with corresponding probabilities p(q,r) and 1-p(q,r). 

 USAGE: ENTR = h_func(P)

 INPUT:

 P
         A QxR matrix of the probabilities p(q,r), q = 1,2,...,Q, 
         r = 1,2,...,R. 

 OUTPUT:

 ENTR    
         A QxR matrix with its (q,r) element being the entropy of 
         the binary random variable X(q,r).
   
"""
function h_func(p)

  entVal = fill(0.0,size(p))
  @inbounds for i in eachindex(p,entVal)
    if p[i]>0 && p[i]<1
      entVal[i] = -1 * p[i]*log2( p[i] ) -
                  (1 - p[i])*log2( 1 -p[i] )
    end
  end
  return entVal
end
