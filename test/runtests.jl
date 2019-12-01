#informME: An information-theoretic pipeline for WGBS data
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

using InformMe
using BlackBoxOptim
using Test
import Statistics.mean
using FileIO
using JLD2

@testset "Test Core Algs" begin
    an = zeros(Float64,10)
    cn = zeros(Float64,9)

    logZ1,logZ0,logZ = InformMe.computeZ(an,cn)
    @test logZ ≈ (10*log(2))

    p1,transProbs = InformMe.computeMCtransProbs(an,cn,logZ1,logZ0,logZ)
    @test p1 ≈ 0.5
    @test transProbs ≈ 0.5.*ones(9,2)

    logZ1tilde,logZ0tilde,logZtilde = InformMe.computeZtilde(an,cn)
    @test logZtilde ≈ (10*log(2))

    logMargProb = InformMe.calcMargProb(1,2,Int8[0,1,0],logZ1,logZ0,logZ,logZ1tilde,logZ0tilde,an,cn)
    @test logMargProb ≈ 3*log(0.5)

    yProbs,yVals,NMEy = InformMe.computeLstats(0.5,0.5*ones(10,2))
    @test yProbs ≈  [0.00048828125,
                     0.00537109375,
                     0.02685546875,
                     0.08056640625,
                     0.1611328125,
                     0.2255859375,
                     0.2255859375,
                     0.1611328125,
                     0.08056640625,
                     0.02685546875,
                     0.00537109375,
                     0.00048828125]
    @test NMEy ≈ 0.7742121153511642

    yProbs,yVals,NMEy = InformMe.computeLstats(0.5,0.5*ones(20,2))
    @test NMEy ≈ 0.72805193277145

    x,y,z=InformMe.maxent([0.5],0:(1/10):1)
    @test y ≈ 0.09090909090909093.*ones(11)

    Xsamp = InformMe.exactSampling(0.5,0.5*ones(1000000,2))
    @test 0.49 <= mean(Xsamp) <= 0.51

    C,ia,ic = InformMe.uniqueOct(["test","testing","test","testing","testtest"])
    @test C == ["test","testing","testtest"]
    @test ia == [1,2,5]
    @test ic == [1,2,1,2,3]

    entVals = InformMe.h_func([-1 1.0;0.5 0.75])
    @test entVals ≈ [0.0 0.0;1.0 0.8112781244591328]

    i,j = InformMe.findSortedIndices(collect(-1:2:100),-5,1)
    @test [i,j] == [1,2]
    i,j = InformMe.findSortedIndices(collect(-1:2:100),100,101)
    @test [i,j] == [0,-1]
    i,j = InformMe.findSortedIndices(collect(-1:2:100),99,101)
    @test [i,j] == [51,51]
    i,j = InformMe.findSortedIndices(collect(-1:2:100),-5,-4)
    @test [i,j] == [0,-1]
    i,j = InformMe.findSortedIndices(collect(-1:2:100),1,5)
    @test [i,j] == [2,4]

    dens = InformMe.nDensity(collect(1:200:401),2,400)
    @test dens ≈ (3.0/400)

    density = ones(Float64,10)
    Dist = 2 .* ones(Int64,9)
    an,cn = InformMe.computeAnCn(ones(10),Dist,[1.;2;3;4;5])
    @test an ≈ [4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0]
    @test cn ≈ [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]

    an=zeros(10);cn=zeros(9);
    aveLogLike=InformMe.computeAveLogLikelihood(an,cn,zeros(Int8,8,10)',
                                                ones(Int,8),10 .*ones(Int,8))
    @test aveLogLike ≈ 10*log(0.5)

    an=ones(10)
    aveLogLike=InformMe.computeAveLogLikelihood(an,cn,ones(Int8,10,10)',
                                                ones(Int,10),10 .*ones(Int,10))
    @test aveLogLike ≈ -1.2692801104297242

    cn=ones(9)
    aveLogLike=InformMe.computeAveLogLikelihood(an,cn,ones(Int8,10,10)',
                                                ones(Int,10),10 .*ones(Int,10))
    @test aveLogLike ≈ -0.06417626486765826

    an=collect(1:10.)
    aveLogLike=InformMe.computeAveLogLikelihood(an,cn,ones(Int8,10,10)',
                                                ones(Int,10),10 .*ones(Int,10))
    @test aveLogLike ≈ -0.018862691685896493

    newMatrix, CpGstart, CpGend = InformMe.processMatrix(Int8[-1 0 -1 1 -1 1;])
    @test newMatrix == [-1   0  -1  -1  -1  -1;
                        -1  -1  -1   1  -1  -1;
                        -1  -1  -1  -1  -1   1]
    @test CpGstart==CpGend==[2,4,6]

end

@testset "Test Larger Functions" begin
    # comptue estimate in a=0,c=0 data
    regionStruct=InformMe.estimateParams(
      "blah",
      "blah",
      10 .* ones(Int64,10),
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      Int8[(rand()>0.5 ? 0 : 1) for ind1=1:100, ind2=1:10])

    # Note model is unidentifiable, so test an,cn not theta
    an,cn = InformMe.computeAnCn(
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      10 .* ones(Int64,10),regionStruct.thetabest)
    @test abs(mean(an)) < 0.1
    @test abs(mean(cn)) < 0.1

    # comptue estimate in a=0,c=inf data
    dataMat=[Int8[1 for ind1=1:100, ind2=1:10];
             Int8[0 for ind1=1:100, ind2=1:10]]
    dataMat[1:5,1:5].= -1
    dataMat[101:105,1:5].= -1
    regionStruct=InformMe.estimateParams(
      "blah",
      "blah",
      10 .* ones(Int64,10),
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      dataMat)

    # Note model is unidentifiable, so test an,cn not theta
    an,cn = InformMe.computeAnCn(
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      10 .* ones(Int64,10),regionStruct.thetabest)
    @test abs(mean(an)) < 0.1
    @test mean(cn) > 1

    # comptue estimate in a=inf,c=inf data
    regionStruct=InformMe.estimateParams(
      "blah",
      "blah",
      10 .* ones(Int64,10),
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      Int8[1 for ind1=1:100, ind2=1:10])

    # Note model is unidentifiable, so test an,cn not theta
    an,cn = InformMe.computeAnCn(
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      10 .* ones(Int64,10),regionStruct.thetabest)
    @test mean(an) > 9
    @test mean(cn) > 1


    # comptue estimate in a=-inf,c=inf data
    regionStruct=InformMe.estimateParams(
      "blah",
      "blah",
      10 .* ones(Int64,10),
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      Int8[0 for ind1=1:100, ind2=1:10])

    # Note model is unidentifiable, so test an,cn not theta
    an,cn = InformMe.computeAnCn(
      Float64[InformMe.nDensity(collect(10:10:100),n,500) for n=1:10],
      10 .* ones(Int64,10),regionStruct.thetabest)
    @test mean(an) < -9
    @test mean(cn) > 1

    # Run fastaToCpG on toy data
    outdir = @__DIR__
    InformMe.fastaToCpG("toy_genome.fasta",outdir=outdir)
    for chr = 1:5
        chr_num_str = string(chr)
        CpGdataT = joinpath(outdir,string("CpGlocationChr", chr_num_str, "_true.jld2"))
        CpGdata = joinpath(outdir,string("CpGlocationChr", chr_num_str, ".jld2"))

        # test each output
        finalCpGloc  = FileIO.load(CpGdata,"finalCpGloc")
        finalCpGlocT = FileIO.load(CpGdataT,"finalCpGloc")
        @test finalCpGloc == finalCpGlocT

        Dist         = FileIO.load(CpGdata,"Dist")
        DistT        = FileIO.load(CpGdataT,"Dist")
        @test Dist == DistT

        density      = FileIO.load(CpGdata,"density")
        densityT     = FileIO.load(CpGdataT,"density")
        @test density ≈ densityT

        CpGlocation  = FileIO.load(CpGdata,"CpGlocation")
        CpGlocationT = FileIO.load(CpGdataT,"CpGlocation")
        @test CpGlocation == CpGlocationT

        # remove test results JLD files
        rm(CpGdata)
    end

    # test matrix from reads
    bamFile = joinpath(@__DIR__,"toy_normal_pe.bam")
    region_str="chr1:1-130"
    command = ["samtools", "view","-q", "30", "-F", "3328", "-f", "3",
                    bamFile, region_str]
    SAMreads = split(readchomp(`$command`),'\n')
    CpGlocation = FileIO.load(joinpath(@__DIR__,"CpGlocationChr1_true.jld2"),"CpGlocation")
    ansMat = InformMe.matrixFromReads(SAMreads,CpGlocation[1:10],true,0)
    @test ansMat == [ -1  -1  -1  -1  -1   1  1  1   1   1;
                      -1   1   0   1   1   1  0  1   1   0;
                      -1   0   1   1   1   1  1  1   1   0;
                      -1  -1  -1  -1  -1  -1  1  1   1   0;
                      -1  -1  -1  -1   0   1  1  1  -1   1;
                      -1  -1   1   1   1   1  1  1   0   1;
                      -1  -1  -1   1   1   1  0  0   1  -1]

    # test parTaskMatrix
    result = InformMe.parTaskMatrix(1,CpGlocation,bamFile,true,0,"1",5,130,true)
    @test result.observedMatrix == ansMat
    
    # test estimation    
    regionStruct=InformMe.estimateParams(
      "blah",
      "blah",
      10 .* ones(Int64,10),
      Float64[0.06 for n=1:10],
      ansMat;boundaryConditions=true)
    an,cn = InformMe.computeAnCn(
      Float64[0.06 for n=1:10],
      10 .* ones(Int64,10),regionStruct.thetabest)

    @test an[2] > 0.5
    @test abs(an[10]) < 0.15

    # test the objective function
    newMatrix, CpGstart, CpGend = InformMe.processMatrix( ansMat )
    @test newMatrix == [ -1  -1  -1  -1  -1   1   1   1   1   1;
                         -1   1   0   1   1   1   0   1   1   0;
                         -1   0   1   1   1   1   1   1   1   0;
                         -1  -1  -1  -1  -1  -1   1   1   1   0;
                         -1  -1  -1  -1   0   1   1   1  -1  -1;
                         -1  -1  -1  -1  -1  -1  -1  -1  -1   1;
                         -1  -1   1   1   1   1   1   1   0   1;
                         -1  -1  -1   1   1   1   0   0   1  -1]    
    @test CpGstart == [6;2;2;7;5;10;3;4]
    @test CpGend == [10;10;10;10;8;10;10;9] 

    DistInRegion,densityInRegion= 10 .* ones(Int64,10), Float64[0.06 for n=1:10]
    thetaBest = [-0.019,13.1767,0.0335,-3.0,-0.0020]
    function objFnToMinimize(theta::AbstractArray)
      An,Cn = InformMe.computeAnCn(densityInRegion,DistInRegion[1:(end-1)],theta)
      aveLogLikelihood::Float64 = InformMe.computeAveLogLikelihood(An,Cn,newMatrix', CpGstart, CpGend)
      return -1*aveLogLikelihood
    end
    obj = objFnToMinimize(thetaBest)
    @test abs(obj-2.8385) < 0.001
end

@testset "Test Full Pipeline" begin
  # clean and setup testing environment
  outdir = @__DIR__

  isdir(joinpath(outdir,"estimation")) && rm(joinpath(outdir,"estimation"),recursive=true)
  isdir(joinpath(outdir,"genome")) && rm(joinpath(outdir,"genome"),recursive=true)
  isdir(joinpath(outdir,"indexedBAMfiles")) && rm(joinpath(outdir,"indexedBAMfiles"),recursive=true)
  isdir(joinpath(outdir,"matrices")) && rm(joinpath(outdir,"matrices"),recursive=true)
  isdir(joinpath(outdir,"output")) && rm(joinpath(outdir,"output"),recursive=true)
  isdir(joinpath(outdir,"scratch")) && rm(joinpath(outdir,"scratch"),recursive=true)

  mkpath(joinpath(outdir,"estimation"))
  mkpath(joinpath(outdir,"genome"))
  mkpath(joinpath(outdir,"indexedBAMfiles"))
  mkpath(joinpath(outdir,"matrices"))
  mkpath(joinpath(outdir,"output"))
  mkpath(joinpath(outdir,"scratch"))

  !islink(joinpath(outdir,"genome","CpGlocationChr1.jld2")) && 
    symlink(joinpath(outdir,"CpGlocationChr1_true.jld2"),joinpath(outdir,"genome","CpGlocationChr1.jld2"))
  !islink(joinpath(outdir,"genome","CpGlocationChr2.jld2")) &&
    symlink(joinpath(outdir,"CpGlocationChr2_true.jld2"),joinpath(outdir,"genome","CpGlocationChr2.jld2"))
  !islink(joinpath(outdir,"genome","CpGlocationChr3.jld2")) &&
    symlink(joinpath(outdir,"CpGlocationChr3_true.jld2"),joinpath(outdir,"genome","CpGlocationChr3.jld2"))
  !islink(joinpath(outdir,"genome","CpGlocationChr4.jld2")) &&
    symlink(joinpath(outdir,"CpGlocationChr4_true.jld2"),joinpath(outdir,"genome","CpGlocationChr4.jld2"))
  !islink(joinpath(outdir,"genome","CpGlocationChr5.jld2")) &&
    symlink(joinpath(outdir,"CpGlocationChr5_true.jld2"),joinpath(outdir,"genome","CpGlocationChr5.jld2"))
  !islink(joinpath(outdir,"indexedBAMfiles","toy_cancer_pe.bam")) && 
    symlink(joinpath(outdir,"toy_cancer_pe.bam"),joinpath(outdir,"indexedBAMfiles","toy_cancer_pe.bam"))
  !islink(joinpath(outdir,"indexedBAMfiles","toy_cancer_pe.bam.bai")) &&
    symlink(joinpath(outdir,"toy_cancer_pe.bam.bai"),joinpath(outdir,"indexedBAMfiles","toy_cancer_pe.bam.bai"))
  !islink(joinpath(outdir,"indexedBAMfiles","toy_normal_pe.bam")) &&
    symlink(joinpath(outdir,"toy_normal_pe.bam"),joinpath(outdir,"indexedBAMfiles","toy_normal_pe.bam"))
  !islink(joinpath(outdir,"indexedBAMfiles","toy_normal_pe.bam.bai")) && 
    symlink(joinpath(outdir,"toy_normal_pe.bam.bai"),joinpath(outdir,"indexedBAMfiles","toy_normal_pe.bam.bai"))
 
  
  # run test
  InformMe.convertBAMtoBits("toy_normal_pe","normal",chr_nums=1:2,
                            numBasesToTrim=10,
                            boundaryConditions=true,
                            reference_path=joinpath(outdir,"genome"),
                            bamfile_path=joinpath(outdir,"indexedBAMfiles"),
                            matrices_path=joinpath(outdir,"matrices"),
                            estimation_path=joinpath(outdir,"estimation"),
                            outdir=joinpath(outdir,"output"))
 
  # test matrices generation
  data = FileIO.load(joinpath(outdir,"matrices","chr1","toy_normal_pe_matrices.jld2"),"mapObjData")
  pbam = data["chr1/bp1-3000"]
  mat  = pbam.observedMatrix
  vect = (sum(mat.>0,dims=1)./sum(mat.>-1,dims=1)) # observed meth freqs
  vect = vect[:]
  obsPerc = mean(vect[.!(isnan.(vect))])
  @test (obsPerc > 0.79) && (obsPerc < 0.81) 

  # test final answer 
  outfile = joinpath(outdir,"output","MML-normal.bed")
  result = parse(Float64,chomp(read(`awk '{if(NR>1){total+=$4}}END{print total/NR}' $outfile`,String)))
  @test (result > 0.77) && (result < 0.83)

end
