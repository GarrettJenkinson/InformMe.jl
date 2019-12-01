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
function analysisParallelTask(startBP, regionSize,chr_num_str,mapObjEstData,CpGlocation,subRegionSize;MSIflag=false,ESIflag=false,MCflag=false)

  #  try

  endBP = startBP+regionSize-1 # the last region might be too large in BP but this does not matter

  # find relative path of this region
  locationPathName = string("chr", chr_num_str, "/bp", startBP, "-", endBP)

  if haskey(mapObjEstData,locationPathName)
    # find CpG sites
    lower_index,upper_index = findSortedIndices(CpGlocation,startBP,endBP)
    CpGlocs_local = CpGlocation[lower_index:upper_index]

    localEstStruct = mapObjEstData[locationPathName]

    # estimate parameters for this region
    FDregionStruct = methAnalysisForRegion(localEstStruct,CpGlocs_local,startBP,endBP,subRegionSize,MSIflag=MSIflag,ESIflag=ESIflag,MCflag=MCflag)
    return FDregionStruct
  else
    return AnalyRegStruct()
  end

  #   catch
  #     println(string("Error in methAnalysisForRegion at ", locationPathName))
  #   end

end
