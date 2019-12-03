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
 This function processes a BAM file with aligned reads to a reference
 genome and produces methylation information for nonoverlapping genomic
 regions (containing the same number of base pairs) in a given chromosome.
 The final output for each genomic region is a matrix with -1,0,1 values.
 Each row of the matrix is a methylation read, whereas each column
 represents a CpG site within the genomic region. A value of -1 indicates
 no methylation information is available for the CPG site, 0 indicates
 that the CpG site is unmethylated, and 1 indicates that the CpG site
 is methylated.

 This function depends on a working instalation of SAMtools that is on
 the system path $PATH.

 Before running this function, FastaToCpG.m must be run ONCE.

## USAGE (default):

 `matrixFromBam(bam_prefix,chr_num)`

## USAGE (optional):

 Example of optional usage with additional input parameters.
 
 `matrixFromBam(bam_prefix,chr_num; reference_path="/path/to/ref")`

## MADATORY INPUTS:

 `bam_prefix`

                Prefix of the BAM file (without the .bam extension). This
                file must be sorted from the least to the greatest base
                pair position along the reference sequence and must be
                indexed (i.e., the associated BAI file must be available).
                The file name must not contain "." characters, but can
                contain "_" instead. Moreover, the file name should be
                unique.

 `chr_num`

                Number representing the chromosome to be processed.

## OPTIONAL INPUTS:

 `reference_path`

                Path to the root subdirectory where the outputs of this
                function are stored.
                Default value: "./genome/"

 `bamfile_path`

                Path to the subdirectory where the BAM file is located.
                Default value: "./indexedBAMfiles/"

 `matrices_path`

                Path to the subdirectory where the output of this function
                is stored.
                Default value: "./matrices/"

 `pairedEnds`

                Flag for paired end read support. A value of true indicates
                that the sequencer employed paired end reads, whereas a
                value of false indicates that the sequencer employed single
                end reads.
                Default value: true

 `numBasesToTrim`

                A vector of integers specifying how many bases should be
                trimmed from the begining of each read. If the vector
                contains two integers, then the first integer specifies
                how many bases to trim from the first read in a read pair,
                whereas the second integer specifies how many bases should
                be trimmed from the second read in the pair. If the
                vector contains one integer, then all reads will have
                that number of bases trimmed from the beginning of the
                read. If no bases are to be trimmed, then this input
                must be set to 0.
                Default value: 0

 `regionSize`

                The size of the genomic regions for which methylation
                information is produced (in number of base pairs).
                Default value: 3000

 `minCpGsReqToModel`

                The minimum number of CpG sites within a genomic region
                required to perform statistical estimation.
                Default value: 10

 The default values of regionSize and minCpGsReqToModel should only be
 changed by an expert with a detailed understanding of the code and the
 methods used.
"""
function matrixFromBam(bamFilename,chr_num;
                       reference_path="./genome/",
                       bamfile_path="./indexedBAMfiles/",
                       matrices_path="./matrices/",
                       pairedEnds=true,
                       numBasesToTrim=0,
                       regionSize=3000,
                       minCpGsReqToModel=10,
                       numProcessors=1)

  # Manual checks/corrections of inputs
  if bamfile_path[end:end]!="/"
    bamfile_path=string(bamfile_path, "/")
  end

  if matrices_path[end:end]!="/"
    matrices_path=string(matrices_path, "/")
  end

  if string(reference_path[end:end])!="/"
    reference_path=string(reference_path, "/")
  end

  ###########################################################################
  # Initialize
  ###########################################################################
  bamFile           = string(bamfile_path, bamFilename, ".bam")
  CpGfileNameRoot   = "CpGlocationChr" # root file name of .mat files
  newLineChar       = '\n'    # constant used for text processing

  #
  # get name of chromsome as a string
  #
  chr_str = string(chr_num)

  #
  # Check if the final output already exists, and exit with error if so
  #

  if isfile(string(matrices_path,"chr", chr_str,
                             "/", bamFilename, "_matrices.jld2"))
    println("Exiting. Final merged output file already exists:")
    println(string(matrices_path,"chr", chr_str,
                               "/", bamFilename, "_matrices.jld2"))
    return
  end


  #
  # load chromosome CpG location information
  #
  CpGlocation=FileIO.load( string(reference_path,
                                  CpGfileNameRoot, chr_str, ".jld2"),
                          "CpGlocation" )
  includeChrInRef=FileIO.load( string(reference_path,
                                      CpGfileNameRoot, chr_str, ".jld2"),
                              "includeChrInRef" )

  numCpGs = length(CpGlocation)

  #
  ###########################################################################
  # Initialize hashtable
  ###########################################################################
  #

  mapObjData = Dict{String,processedBAM}()

  #
  ###########################################################################
  # Proceed through all regions in a chromosome in parallel
  ###########################################################################
  #
  finalCpGloc=CpGlocation[numCpGs]
  allStartBPs = 1:regionSize:finalCpGloc

  # AllData = Array{EstRegStruct}(undef,length(allStartBPs))
  # for x = 1:length(allStartBPs)
  #  AllData[x] = parTaskMatrix(allStartBPs[x],CpGlocation,bamFile,pairedEnds,
  #                             numBasesToTrim,chr_str,minCpGsReqToModel,
  #                             regionSize,includeChrInRef)
  # end
  numProcessors > 1 && addprocs(numProcessors-1)
  numProcessors > 1 && @everywhere eval(quote using informme end)
  AllData = pmap( x -> parTaskMatrix(x,
                              CpGlocation,bamFile,pairedEnds,
                              numBasesToTrim,chr_str,minCpGsReqToModel,
                              regionSize,includeChrInRef),
                allStartBPs )
  numProcessors > 1 && rmprocs(workers())


  # save output to hashtable with correct keys (and only takes values when observedMatrix is present)
  for index = 1:length(allStartBPs)
    if !isempty(AllData[index])
      keyStr = string("chr", chr_str, "/bp", allStartBPs[index], "-", (allStartBPs[index]+regionSize-1))
      mapObjData[keyStr]=AllData[index]
    end
  end

  #
  ###########################################################################
  # Save output to file
  ###########################################################################
  #

  # check that output directory exists, create if not
  if !isdir(string(matrices_path, "chr", chr_str))
    mkpath(string(matrices_path, "chr", chr_str))
  end

  # write out to file
  FileIO.save(string(matrices_path, "chr", chr_str,
                             "/", bamFilename, "_matrices.jld2"),"mapObjData", mapObjData)

  return
end
