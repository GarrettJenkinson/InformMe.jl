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
 Runs the entire MethToBits pipeline (assuming single BAM file as input).
"""
function convertBAMtoBits(bamFilenames,phenoName;
                          reference_path="./genome/",
                          bamfile_path="./indexedBAMfiles/",
                          matrices_path="./matrices/",
                          estimation_path="./estimation/",
                          outdir="./output/",
                          pairedEnds=true,
                          numBasesToTrim=0,
                          minCpGsReqToModel=10,
                          regionSize=3000,
                          boundaryConditions=false,
                          MSIflag=false,
                          ESIflag=false,
                          MCflag=false,
                          subRegionSize=150,
                          chr_nums=1:22,
                          numProcessors=nprocs())


for chr_num in chr_nums
    bamFileNameSplit = split(bamFilenames,",")
    for bamFilename = string.(bamFileNameSplit)
        println(string("Generating matrix from ", bamFilename, " on chromosome ", chr_num))
        matrixFromBam(bamFilename,chr_num,
                       reference_path=reference_path,
                       bamfile_path=bamfile_path,
                       matrices_path=matrices_path,
                       pairedEnds=pairedEnds,
                       numBasesToTrim=numBasesToTrim,
                       regionSize=regionSize,
                       minCpGsReqToModel=minCpGsReqToModel,
                       numProcessors=numProcessors)
    end

    println(string("Estimating model for ", phenoName, " on chromosome ", chr_num))
    estParamsForChr(bamFilenames,phenoName,matrices_path,reference_path,chr_num,
                    estimation_path=estimation_path,
                    regionSize=regionSize,
                    boundaryConditions=boundaryConditions,
                    numProcessors=numProcessors)

    println(string("Analyzing model for ", phenoName, " on chromosome ", chr_num))
    methAnalysisForChr(phenoName,chr_num,reference_path,estimation_path,
                       outdir=outdir,
                       MSIflag=MSIflag,
                       ESIflag=ESIflag,
                       MCflag=MCflag,
                       regionSize=regionSize,
                       subRegionSize=subRegionSize,
                       numProcessors=numProcessors)

 end


  println(string("Making BED files for ", phenoName))
  makeBedsForMethAnalysis(phenoName,outdir,reference_path,
                              outdir=outdir,
                              chrs=[string("chr",i) for i=chr_nums],
                              MSIflag=MSIflag,
                              ESIflag=ESIflag,
                              MCflag=MCflag,
                              regionSize=regionSize,
                              subregionSize=subRegionSize)

  return
end
