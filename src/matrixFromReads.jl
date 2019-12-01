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
 This function determines the methylation status of the CpG sites within
 a given genomic region using SAM-formatted reads with base pair coverage
 within the region. The output is a matrix with -1,0,1 values. Each row
 of the matrix is a single sequencing read, whereas each column represents
 a CpG site within the genomic region. A value of -1 indicates no
 methylation information is available for the CPG site, 0 indicates that
 the CpG site is unmethylated, and 1 indicates that the CpG site is
 methylated.

 USAGE:

 observedMatrix = matrixFromReads(currentReads,CpGsInRegion,...
                                       pairedEnds,numBasesToTrim)

 INPUTS:

 currentReads
                  Cell array containing a list of strings that are
                  SAM-formatted reads with base pair coverage within
                  the present genomic region.

 CpGsInRegion
                  Vector containing the list of the base pair locations
                  of the CpG sites within the present region.

 pairedEnds
                  Flag for paired end read support. A value of 1 indicates
                  that the sequencer employed paired end reads, whereas a
                  value of 0 indicates that the seqeuncer employed single
                  end reads.

 OUTPUT:

 observedMatrix
                  Matrix whose (m,n) element takes one of the following
                  three values:
                   -1 if the m-th read does not observe the n-th CpG site.
                    0 if the m-th read indicates that the n-th CpG site
                      is unmethylated.
                    1 if the m-th read indicates that the n-th CpG site
                      is methylated.
"""
function matrixFromReads(currentReads,CpGsInRegion,pairedEnds,numBasesToTrim)

  numRawReads    = length(currentReads)
  N              = length(CpGsInRegion)
  observedMatrix = -1*ones(Int8,numRawReads,N)
  readNames      = Array{String}(undef,numRawReads)

  #
  # convert reads into (-1,0,1) vectors
  #

  for read_num=1:numRawReads

    # seperate read into its \t (tab) delimited elements
    processedRead = split(currentReads[read_num],"\t")

    # store for later
    readNames[read_num] = processedRead[1]

    #
    # find start and end of read
    #
    startOfRead  = parse(Int,processedRead[4])       # 4-th column is start (leftmost base pair) of read
    lengthOfRead = length(processedRead[10])   # 10-th column is the string of base pairs
    endOfRead    = startOfRead+lengthOfRead-1  # last observed base pair in read

    #
    # only deal with reads that are free of deletions, insertions,
    # and clippings
    #
    if !occursin(r"[idnshp]"i,processedRead[6])
      # make sure only base pair matches or mismatches (no indels, etc)

      # find which strand the read is on (forward 0 or reverse 1)
      SAMflag = bitstring(parse(UInt16,processedRead[2]))[end:-1:5] # 2-nd column is bitwise flag
      # matlab command:   SAMflag = bitget(uint16(str2double(processedRead[2])),1:12,'uint16')

      strand  = SAMflag[5]
      # 5-th flag tells whether the read is reverse complemented


      # find all locations in string that have CpG information

      # nOfFirstCpGinfo is the number (from 1 to N) of first CpG site observed
      # nOfLastCpGinfo is the number of the last CpG site observed on read
      nOfFirstCpGinfo,nOfLastCpGinfo = findSortedIndices(CpGsInRegion,startOfRead,endOfRead-1)

      indicesOfCpGinfoInStr = CpGsInRegion .- startOfRead .+ 1
      # index within the string of methylation information
      # (this can go negative or become greater than the
      # string length for CpG sites not observed)

      # process these CpG sites
      for n = nOfFirstCpGinfo:nOfLastCpGinfo   # indexing of CpG sites

        indexOfInfo = indicesOfCpGinfoInStr[n] # index of CpG site information on string

        informativeBases = processedRead[10][indexOfInfo:(indexOfInfo+1)]
        # read at "CG" position

        # ASCII conversion to integer, subtract 33 to convert to Phred quality
        informativeBaseQualityC = Int(processedRead[11][indexOfInfo])-33
        informativeBaseQualityG = Int(processedRead[11][indexOfInfo+1])-33

        # Check if this location should be trimmed
        if '1'==SAMflag[7] && '0'==SAMflag[8] # first read in a paired end read
          currNumBasesToTrim = numBasesToTrim[1]
        elseif '0'==SAMflag[7] && '1'==SAMflag[8] # second read in a paired end read
          currNumBasesToTrim = numBasesToTrim[end]
        else # not a paired end read
          currNumBasesToTrim = numBasesToTrim[1]
        end

        if '0'==strand && (indexOfInfo<=currNumBasesToTrim) # forward strand and beginning of mapped read
          continue # trim this location at front of mapped read
        elseif '1'==strand && ((indexOfInfo+1)>=(lengthOfRead-currNumBasesToTrim+1)) # reverse strand and end of mapped read
          continue # trim this location at back of mapped read
        end

        # Gather information from location if high quality
        if (informativeBaseQualityC >= 20) && (informativeBaseQualityG >= 20)
          # this is a reliable piece of
          # information (1 in 100 chance
          # of error) and it is not in
          # the first 1:numBasesToTrim
          # basepairs of the read

          #
          # find methylation information (if present)
          #

          if occursin(r"CG"i,informativeBases)#strcmpi(informativeBases,"CG")
            observedMatrix[read_num,n] = 1 # methylation
          elseif occursin(r"TG"i,informativeBases)||occursin(r"CA"i,informativeBases)
            observedMatrix[read_num,n] = 0 # no methylation
          end # end of comparing informativeBases

        end # end of processing reliable base pair

      end # end of processing CpG sites from this read

    end # end of perfect alignment match

  end # end of processing all reads

  #
  # find which reads actually contain CpG information
  #
  observationUseful = vec(sum(observedMatrix.>-1,dims=2).>0)

  #
  # create final matrix (merging paired ends if required)
  #
  if pairedEnds

    observedMatrixUsefulSingleReads = observedMatrix[observationUseful,:] # filter out non-useful observations

    readNamesUsefulSingleReads=readNames[observationUseful]

    # find non-unique readnames (i.e., read pairs)
    uniqueReads,IA,IC = uniqueOct(readNamesUsefulSingleReads) # find read pairs

    observedMatrix = observedMatrixUsefulSingleReads[IA,:] # all "first" pairs

    for index=1:length(IC)
      for col=1:length(observedMatrix[IC[index],:])
        observedMatrix[IC[index],col] = max(observedMatrix[IC[index],col],
                                            observedMatrixUsefulSingleReads[index,col])
      end
                                    # clever use of max function:
                                    # -1 replaced by 0 or 1 (observation replaces unobserved)
                                    #  0 replaced by 1 (hemimethylated is methylated)
                                    # note repeated use of max on same row has no effect

    end

  else # no paired ends to combine, so just remove reads that did
       # not provide methylation state
    observedMatrix = observedMatrix[observationUseful,:]
  end

  return observedMatrix

end
