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
 This function takes a list of BAM files (which correspond to the same
 phenotype) and performs statistical model estimation within a specific
 chromosome of interest. The function can be used on a computing cluster
 to break the work of model estimation to many independent parallel job
 processes. This is performed only after matrixFromBam.m.

 USAGE (default):

 estParamsForChr(mat_files,prefix,matrices_path,reference_path,chr_num)

 USAGE (optional):

 Example of optional usage with additional input parameters.
 estParamsForChr(mat_files,prefix,matrices_path,reference_path,chr_num,
                 regionSize=2000)

 MANDATORY INPUTS:

 mat_files
        All the .mat files to be included in the model. This can be a
        single .mat file or multiple files in the form of a comma-sepa-
        rated list of files.

 prefix
        A string that specifies the name of the modeled phenotype.
        The output files produced will contain this prefix.

 matrices_path
        A string that specifies the path to the directory that
        where the output will be stored.

 reference_path
        A string that specifies the path to the directory that
        contains the results of analysis of the reference genome
        performed by FastaToCpG.m as well as the results of
        methylation calling performed by matrixFromBam.m.

 chr_num
        Chromosome number 1 to 22 (in humans) specifying the
        chromosome for which statistical estimation must be
        performed.

 OPTIONAL INPUTS:

 regionSize
        The size of the genomic region used for parameter
        estimation (in number of base pairs).
        Default value: 3000

 boundaryConditions
        Flag to decide if boundary conditions should be estimated
        freely in MLE.
        Default value: false

 The default value of regionSize should only be changed by an expert with
 a detailed understanding of the code and the methods used.
"""
function estParamsForChr(mat_files,phenoName,matrices_path,reference_path,chr_num;
                         estimation_path="./results/",
                         regionSize=3000,
                         boundaryConditions=false,
                         numProcessors=1)

  # Manual checks/corrections of inputs
  if matrices_path[end:end]!="/"
    matrices_path=string(matrices_path, "/")
  end
  if reference_path[end:end]!="/"
    reference_path=string(reference_path,"/")
  end
  if estimation_path[end:end]!="/"
    estimation_path=string(estimation_path,"/")
  end

  chr_num_str = string(chr_num)

  ###########################################################################
  # Check if the final output already exists, and exit with error if so
  ###########################################################################
  outfile = string(estimation_path, "chr", chr_num_str, "/",phenoName,"_fit.jld2")
  if isfile(outfile)
    println("Exiting. Final output file already exists:")
    println(outfile)
    return
  end

  #
  ###########################################################################
  # Find which regions are assigned to this processor
  ###########################################################################
  #

  # Load CpG site data on chromosome (to determine how many regions to loop through, etc.)
  CpGdata = string(reference_path, "CpGlocationChr", chr_num_str, ".jld2")
  finalCpGloc = FileIO.load(CpGdata,"finalCpGloc")
  Dist        = FileIO.load(CpGdata,"Dist")
  density     = FileIO.load(CpGdata,"density")
  CpGlocation = FileIO.load(CpGdata,"CpGlocation")

  # find all start base-pairs for regions to be modeled on this chromsome
  allStartBPs  = 1:regionSize:finalCpGloc # only need to go to the last CpG site, not last BP

  #
  ###########################################################################
  # Load data for this chromosome
  ###########################################################################
  #
  bamFileNames = split(mat_files,",")
  dataMapObj = Array{Dict{String,processedBAM}}(undef,length(bamFileNames))
  datFileCount = 0
  for dataFile = String.(bamFileNames)
    datFileCount += 1
    fileName = string(matrices_path,"chr",chr_num_str,"/",
                      dataFile,"_matrices.jld2")
    dataMapObj[datFileCount] = FileIO.load(fileName,"mapObjData")
  end

  #
  ###########################################################################
  # Loop through all regions on chromosome in parallel
  ###########################################################################
  #

  # AllData = Array{EstRegStruct}(undef,length(allStartBPs))
  # for x = 1:length(allStartBPs)
  #  AllData[x] = estParamsParallelTask(allStartBPs[x],regionSize,dataMapObj,CpGlocation,Dist,density,phenoName,chr_num_str,boundaryConditions)
  # end
  numProcessors > 1 && addprocs(numProcessors-1)
  numProcessors > 1 && @everywhere eval(quote using informme end)
  AllData = pmap(x -> estParamsParallelTask(x,
                      regionSize,dataMapObj,CpGlocation,Dist,density,phenoName,
                      chr_num_str,boundaryConditions),
                 allStartBPs)
  numProcessors > 1 && rmprocs(workers())


  #initialize hashtable
  mapObjEstData = Dict{String,EstRegStruct}()

  # save output to hashtable with correct keys (and only takes values when observedMatrix is present)
  for index=1:length(allStartBPs)
    if !isempty(AllData[index])
      keyStr = string("chr", chr_num_str, "/bp", allStartBPs[index], "-", (allStartBPs[index]+regionSize-1))
      mapObjEstData[keyStr]=AllData[index]
    end
  end

  # Save output to file

  # check that output directory exists, create if not
  if !isdir(string(estimation_path, "chr", chr_num_str))
    mkpath(string(estimation_path, "chr", chr_num_str))
  end

  # write out to file
  FileIO.save(string(estimation_path, "chr", chr_num_str, "/", phenoName, "_fit.jld2"),
           "mapObjEstData", mapObjEstData)
  return
end
