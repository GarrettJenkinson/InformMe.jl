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
@doc raw"""
 This function performs methylation analysis of a given chromosome in a
 single phenotype. The function can be used on a computing cluster to
 break the analysis work to many independent parallel job processes.
 This is performed only after estParamsForChr.m in the Modeling
 subdirectory is run to build the Ising models for the phenotype.

 USAGE (default):

 methAnalysisForChr(prefix,chr_num,reference_path,estimation_path)

 USAGE (optional):

 Example of optional usage with additional input parameters.
 methAnalysisForChr(prefix,chr_num,reference_path,estimation_path,
                    outdir="/path/to/output")

 MANDATORY INPUTS:

 prefix
               A string that specifies the name of the phenotype to be
               analyzed.

 chr_num
               Chromosome number (1 to 22 in humans) specifying the
               chromosome for which methylation analyis must be
               performed.

 reference_path
               A string that specifies the path to the directory that
               contains the results of analysis of the reference genome
               performed by FastaToCpG.m as well as the results of
               methylation calling performed by matrixFromBam.jl.

 estimation_path
               A string that specifies the path to the directory that
               contains the results of parameter estimation performed
               by estParamsForChr.jl.

 OPTIONAL INPUTS:

 outdir
               A string that specifies the path of the directory in which
               the methylation analysis results are written.
               Default value: "./results/"

 MSIflag
               Flag that determines whether this function performs
               computation of the methylation sensitivity index (MSI).
               false: no MSI computation.
               true: allow MSI computation.
               Default value: false

 ESIflag
               Flag that determines whether this function performs
               computation of the entropic sensitivity index (ESI).
               false: no ESI computation.
               true: allow ESI computation.
               Default value: false

 MCflag
               Flag that determines whether this function performs
               computation of turnover ratios, CpG entropies, capacities,
               and relative dissipated energies of methylation
               channels (MCs).
               false: no MC computations.
               true: allow MC computations.
               Default value: false

 regionSize
               The size of the genomic regions used for parameter
               estimation (in number of base pairs).
               Default value: 3000

 subregionSize
               The size of the subregions of a genomic region used
               for methylation analysis (in number of base pairs).
               The ratio regionSize/subregionSize must be an integer.
               Default value: 150

 The default values of regionSize and subregionSize should only be
 changed by an expert with a detailed understanding of the code and
 the methods used.

"""
function methAnalysisForChr(phenoName,chr_num,reference_path,estimation_path;
                        outdir="./results/",
                        MSIflag=false,
                        ESIflag=false,
                        MCflag=false,
                        regionSize=3000,
                        subRegionSize=150,
                        numProcessors=1)


  #
  # Manual checks/corrections of inputs
  #

  if string(reference_path[end:end])!="/"
    reference_path=string(reference_path, "/")
  end
  if string(estimation_path[end:end])!="/"
    estimation_path=string(estimation_path, "/")
  end
  if string(outdir[end:end])!="/"
    outdir=string(outdir, "/")
  end

  # convert chr_num to string
  chr_num_str=string(chr_num)

  #
  ###########################################################################
  # Check if the final output already exists, and exit with error if so
  ###########################################################################
  #
  merged_file=string(outdir,"chr",chr_num_str,"/",phenoName,"_analysis.jld2")
  if isfile(merged_file)
    println("Exiting. Final output file already exists:")
    println(merged_file)
    return
  end

  ###########################################################################
  # load CpG data
  ###########################################################################
  CpGdata     = string(reference_path,"/", "CpGlocationChr", chr_num_str, ".jld2")
  finalCpGloc = FileIO.load(CpGdata,"finalCpGloc")
  Dist        = FileIO.load(CpGdata,"Dist")
  density     = FileIO.load(CpGdata,"density")
  CpGlocation = FileIO.load(CpGdata,"CpGlocation")


  # Find which regions are assigned to this processor
  # find all start base-pairs for regions to be modeled on this chromsome
  allStartBPs = 1:regionSize:finalCpGloc # only need to go to the last CpG site, not last BP

  #
  ###########################################################################
  # Load estimation data for this chromosome
  ###########################################################################
  #

  EstData       = string(estimation_path, "chr", chr_num_str,
                         "/", phenoName, "_fit.jld2")
  mapObjEstData = FileIO.load(EstData,"mapObjEstData")

  #
  ###########################################################################
  # loop over regions in parallel
  ###########################################################################
  #
  # AllData = Array{EstRegStruct}(undef,length(allStartBPs))
  # for x =1:length(allStartBPs)
  #  AllData[x] = analysisParallelTask(allStartBPs[x], regionSize,chr_num_str,mapObjEstData,CpGlocation,subRegionSize,MSIflag=MSIflag,ESIflag=ESIflag,MCflag=MCflag)
  # end
  numProcessors > 1 && addprocs(numProcessors-1)
  numProcessors > 1 && @everywhere eval(quote using informme end)
  AllData = pmap( x -> analysisParallelTask(x,
                          regionSize,chr_num_str,mapObjEstData,CpGlocation,
                          subRegionSize,MSIflag=MSIflag,ESIflag=ESIflag,
                          MCflag=MCflag),
                  allStartBPs)
  numProcessors > 1 && rmprocs(workers())



  #initialize hashtable
  mapObjAnaly = Dict{String,AnalyRegStruct}()

  # save output to hashtable with correct keys (and only takes values when observedMatrix is present)
  for index=1:length(allStartBPs)
    if !isempty(AllData[index])
      keyStr = string("chr", chr_num_str, "/bp", allStartBPs[index], "-", (allStartBPs[index]+regionSize-1))
      mapObjAnaly[keyStr]=AllData[index]
    end
  end

  #
  ###########################################################################
  # Save output to file
  ###########################################################################
  #

  # check that results folder exists, create if not
  if !ispath(outdir)
    mkpath(outdir)
  end
  if !ispath(string(outdir,"chr",chr_num_str))
    mkpath(string(outdir,"chr",chr_num_str))
  end

  # write out to file
  FileIO.save(merged_file,"mapObjAnaly", mapObjAnaly)

 return
end
