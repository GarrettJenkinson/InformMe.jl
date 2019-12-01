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
function estParamsParallelTask(startBP,regionSize,dataMapObj,CpGlocation,Dist,density,phenoName,chr_num_str,boundaryConditions)

  endBP = startBP+regionSize-1 # the last region might be too large in BP but this does not matter

  # find relative path of this region
  locationPathName = string("chr", chr_num_str, "/bp", startBP, "-", endBP)

  # initialize vars
  DistInRegion=Int64[]
  densityInRegion=Float64[]
  dataMat=Array{Int8}(undef, 0, 0)

  #
  # load data for region
  #
  firstObj=true
  for dataFile = 1:length(dataMapObj)
    if haskey(dataMapObj[dataFile],locationPathName)
      regionDataStruct::processedBAM  = dataMapObj[dataFile][locationPathName]
      if firstObj #get values of Dist and density in region
        lower_index,upper_index = findSortedIndices(CpGlocation,startBP,endBP)
        DistInRegion            = Dist[lower_index:upper_index]
        densityInRegion         = density[lower_index:upper_index]
        dataMat::Array{Int8,2}  = regionDataStruct.observedMatrix
        firstObj = false
      else
        dataMat = vcat(dataMat,regionDataStruct.observedMatrix)
      end
    end
  end

  #
  # Estimate parameters
  #
  if !firstObj
    # estimate parameters for this region
    return estimateParams(locationPathName,phenoName,DistInRegion,densityInRegion,dataMat,boundaryConditions=boundaryConditions)
  else
    return EstRegStruct()
  end

end
