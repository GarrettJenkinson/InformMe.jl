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
 This function computes the density of the n-th CpG site within a given 
 chromosome. 

 USAGE:

 density = nDensity(CpGlocation,n,wsize)

 INPUTS: 

 CpGlocation
             Sorted vector of CpG locations within a given chromosome.

 n
             Index of the CpG site within the chromosome.

 wsize
             Window size used in density calculation.

 OUTPUT:

 density
             Computed density of CpG site.

"""
function nDensity(CpGlocation::Array{<:Integer,1},n::Integer,L::Integer)


  density::Float64 = 0.0 # use double precision to represent density

  #
  # compute upper and lower bounds of CpG sites which are included in
  # density calculation
  #
  UpperBound = CpGlocation[n]+floor(L/2)
  LowerBound = CpGlocation[n]-floor(L/2)

  #
  # find the smallest and largest indices of the CpG sites between the
  # upper and lower boundary (inclusive the boundary)
  #
  lower_index,upper_index = findSortedIndices(CpGlocation,LowerBound,UpperBound)
  # CpGlocation must be sorted low to high

  if lower_index==0 # should not happen since current CpG always inside window
    NumberOfCpGsBetweenBounds = 1 # current CpG always inside window
    println("Warning: No CpGs found in density window.")
  else
    NumberOfCpGsBetweenBounds = 1 + upper_index - lower_index # +1 because inclusive at boundaries
  end

  density = NumberOfCpGsBetweenBounds/L
  return density

end
