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

function parTaskMatrix(startBP,CpGlocation,bamFile,pairedEnds,numBasesToTrim,
                       chr_str,minCpGsReqToModel,regionSize,includeChrInRef)

# try


  #                                only need to go to the last CpG site,
  #                                not last base pair
  endBP = startBP+regionSize-1
  #                                last region might be too large in
  #                                base pairs but this does not matter

  #
  # find all CpG sites in the region
  #
  lower_CpGindex,upper_CpGindex = findSortedIndices(CpGlocation,startBP,endBP)
  CpGlocInRegion = CpGlocation[lower_CpGindex:upper_CpGindex]

  #
  # check if there is a sufficient number of CpG sites in region
  #
  if length(CpGlocInRegion) < minCpGsReqToModel
    # skip region since there are not enough CpG sites to model
    return processedBAM()
  end

  #
  #######################################################################
  # Gather all reads relevant to the current region using samtools
  #######################################################################
  #
  # Define the (possibly smaller than 3000 bp) subregion where reads
  # will be collected. Reads with overlap to this region may contain
  # observations of CpG sites within the current 3000 bp genomic region.
  # We only care about the region having CpG sites (this is why not using
  # startBP & endBP). +1 is assigned to the last CpG site so we do not
  # miss a read on the reverse complementary strand.
  #
  if includeChrInRef
    region_str = string("chr", chr_str, ":", CpGlocInRegion[1], "-",(CpGlocInRegion[end]+1))
  else
    region_str = string(chr_str, ":", CpGlocInRegion[1], "-",(CpGlocInRegion[end]+1))
  end

  #
  # Make the samtools view command with the following options
  # -f 3    makes sure only reads that have both paired ends mapped
  #         are counted [3 = (2^0)+(2^1)]
  # -q 30   keeps reads with Phred mapping alignment score >= 30
  #         (i.e., error probability <= 1/1000)
  # -F 3328 excludes PCR duplicates, secondary alignments, and
  #         chimeric alignments [3328 = (2^8)+(2^10)+(2^11)]
  #
  if pairedEnds
    #command      = string("samtools view -q 30 -F 3328 -f 3 ", bamFile, " ", region_str);
    #commandCount = string("samtools view -c -q 30 -F 3328 -f 3 ", bamFile, " ", region_str);
    command      = ["samtools", "view","-q", "30", "-F", "3328", "-f", "3",  bamFile, region_str]
    commandCount = ["samtools", "view","-c","-q", "30", "-F", "3328", "-f", "3",  bamFile, region_str]
  else
    #command      = string("samtools view -q 30 -F 3329 ", bamFile, " ", region_str);
    #commandCount = string("samtools view -q 30 -F 3329 ", bamFile, " ", region_str);
    command      = ["samtools", "view","-q", "30", "-F", "3329",  bamFile, region_str]
    commandCount = ["samtools", "view","-c","-q", "30", "-F", "3329", bamFile, region_str]
    # ignore -f since no paired ends
    # add (2^0) to -F to exclude a paired end read
  end

  #
  # Before running command, check that region is not a overly highly
  # mapped region. If too many reads are mapped to this region, then this
  # indicates that the mapping is unreliable (e.g., due to repeating
  # elements that cannot be uniquely mapped to a reference genome)
  #

  countedNumReadsStr = split(readchomp(`$commandCount`),'\n')

  countedNumReads=parse(Int64,countedNumReadsStr[1])

  if countedNumReads<=0
    return processedBAM()
  elseif countedNumReads>5000
    println(string("Too many reads mapped to region: ", region_str))
    return processedBAM()
  end

  #
  # Now run actual command, given that the number of reads is normal
  #
  SAMreads = split(readchomp(`$command`),'\n')

  if isempty(SAMreads)
    return processedBAM()
  end

  #
  #######################################################################
  # Process the reads relevant to the current region
  #######################################################################
  #

  observedMatrix = matrixFromReads(SAMreads,CpGlocInRegion,pairedEnds,numBasesToTrim)

  #
  #######################################################################
  # Write output
  #######################################################################
  #
  return processedBAM(observedMatrix)



#catch
#   println(string("Error in matrixFromBAM in region start bp", startBP))
#   error()
# end

end
