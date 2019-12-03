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
@doc raw"""
 Runs the entire MethToBits pipeline.

 This function depends on a working instalation of SAMtools that is on
 the system path $PATH.

 Before running this function, FastaToCpG.m must be run ONCE.

## USAGE (default):

 `convertBAMtoBits(bamFilenames,phenoName)`

## USAGE (optional):

 Example of optional usage with additional input parameters.

 `matrixFromBam(bam_prefix,chr_num,reference_path="/path/to/ref")`

## MADATORY INPUTS:

 `bamFilenames`

                A comma seperated string with list of input bam file
                names without the ".bam" extension. These
                files must be sorted from the least to the greatest base
                pair position along the reference sequence and must be
                indexed (i.e., the associated BAI file must be available).
                The file name must not contain "." characters, but can
                contain "_" instead. Moreover, the file name should be
                unique.

 `phenoName`

                A string which will be the unique identifier of this 
                sample/model that is built.


##OPTIONAL INPUTS:

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

 `estimation_path`

                A string that specifies the path to the directory that
                contains the results of parameter estimation performed
                by estParamsForChr.jl.
                Default value: "./estimation/"


 `outdir`

               A string that specifies the path of the directory in which
               the methylation analysis results are written.
               Default value: "./output/"

 `pairedEnds`
         
                Flag for paired end read support. A value of 1 indicates
                that the sequencer employed paired end reads, whereas a
                value of 0 indicates that the sequencer employed single
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

 `boundaryConditions`

                Flag to decide if boundary conditions should be estimated
                freely in MLE.
                Default value: false

 `MSIflag`

               Flag that determines whether this function performs
               computation of the methylation sensitivity index (MSI).
               false: no MSI computation.
               true: allow MSI computation.
               Default value: false

 `ESIflag`

               Flag that determines whether this function performs
               computation of the entropic sensitivity index (ESI).
               false: no ESI computation.
               true: allow ESI computation.
               Default value: false

 `MCflag`

               Flag that determines whether this function performs
               computation of turnover ratios, CpG entropies, capacities,
               and relative dissipated energies of methylation
               channels (MCs).
               false: no MC computations.
               true: allow MC computations.
               Default value: false

 `chr_nums`

              A vector with the chromosomes to be processed (without 
              "chr" string). 
              Default value: 1:22

 `numProcessors`

              The number of processors to use in the computations.
              Note that julia must be started as "julia -p 4" if
              four processors are desired. The nprocs() function
              tells how many cores are available in julia, and
              we default to use them all.
              Default value: nprocs() 

 The default values of regionSize and minCpGsReqToModel should only be
 changed by an expert with a detailed understanding of the code and the
 methods used.

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
