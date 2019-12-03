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
 This function is used to analyze a reference genome in order to find
 and store the locations of all CpG sites within each chromosome and
 to compute the CpG densities at each CpG site as well as the distances
 between neighboring CpG sites. A 1-based coordinate system is used,
 in which the first base is assigned to position 1 and the location
 of a CpG site is defined by the position of the C nucleotide on the
 forward strand of the reference genome.

 This function must be run ONLY ONCE before proceeding with analysis
 of BAM files.

## USAGE (default):

 `fastaToCpg(FASTAfilename)`

## USAGE (optional):

 Example of optional usage with additional input parameters.

 `FastaToCpG(FASTAfilename,outdir="/path/to/outputdir/")`

## MADATORY INPUT:

 `FASTAfilename`

               Full path of FASTA-formatted reference genome to which
               available BAM files have been aligned to.

## OPTIONAL INPUTS:

 `outdir`

               Path where the output will be stored at.
               Default value: "./"

 `wsize`

               Window size used in CpG density calculations.
               Default value: 1000

"""
function fastaToCpG(FASTAfilename;
                    outdir="./",
                    wsize=1000)


  #
  # Manual checks/corrections of inputs
  #

  if outdir[end:end]!="/"
    outdir=string(outdir, "/")
  end

  referenceFASTAgenome = string(outdir, FASTAfilename)

  #
  ###########################################################################
  # nucleotides. Following is not needed: A_int = nt2int('A');
  ###########################################################################
  #
  C_int = BioSequences.DNA_C
  G_int = BioSequences.DNA_G
  T_int = BioSequences.DNA_T

  #
  ###########################################################################
  # sequentially proceed through each chromosome in genome
  ###########################################################################
  #
  fr = BioSequences.FASTA.Reader(open(referenceFASTAgenome, "r"))
  for record in fr
    chr_str = BioSequences.FASTA.identifier(record)

    # convert to string to represent appropriate chromosome
    if (length(chr_str)>3) && (chr_str[1:3]=="chr")
      chr_num_str = chr_str[4:end]
      includeChrInRef = true
    else
      chr_num_str = chr_str
      chr_str = string("chr",chr_num_str)
      includeChrInRef = false
    end

    #
    #######################################################################
    # read seqeuence of current chromosome
    #######################################################################
    #
    sequenceStr = BioSequences.FASTA.sequence(record)
    #
    #######################################################################
    # find all CpG sites for the current chromosome
    #######################################################################
    #
    prevNC  = T_int    # this variable indicates the location of the
    # previous nucleotide - initialized not to
    # be a C
    CpGcount::Int64 = 0
    CpGlocation     = zeros(Int64,Int64(1e7)) # initialize locations of CpG sites
    chrLength       = length(sequenceStr) # compute length of chromosome

    for location = 1:chrLength

      # get location of current nucleotide
      currNC = sequenceStr[location]

      # check for CpG site
      if prevNC == C_int && currNC== G_int
        # CpG found
        # previous location was a CpG site (we consider
        # the C to be the location of the CpG site)
        CpGcount=CpGcount+1
        CpGlocation[CpGcount] = location-1
        # subtract 1 since we are at the G and want location of C
      end
      prevNC=currNC
    end
    CpGlocation = CpGlocation[1:CpGcount]

    #
    #######################################################################
    # make Dist variable for distance
    #######################################################################
    #
    Dist = cat((CpGlocation[2:CpGcount]-CpGlocation[1:(CpGcount-1)]),
               4294967295,dims=1)
    # max value of uint32 ... just chosen as
    # a large number since the distance
    # to the next CpG site is undefined for
    # the last CpG site on a chromosome
    #
    #######################################################################
    # compute density variable
    #######################################################################
    #
    density = [nDensity(CpGlocation,cpgNum,wsize) for cpgNum=1:CpGcount]

    if CpGcount>0
      finalCpGloc = CpGlocation[CpGcount]
    else
      finalCpGloc = NaN
    end
    #
    #######################################################################
    # write output file
    #######################################################################
    #
    FileIO.save(string(outdir, "CpGlocationChr", chr_num_str, ".jld2"),
             "CpGlocation",CpGlocation,"density",density,"Dist",Dist,
             "finalCpGloc",finalCpGloc,"chrLength",chrLength,
             "includeChrInRef",includeChrInRef)
#    println("Chr$chr_num_str.  FinalCpGloc: $finalCpGloc; chrLength: $chrLength; includeChrInRef: $includeChrInRef")
  end

  close(fr)
  return 0
end
