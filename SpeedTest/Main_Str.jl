# cd("Julia/MPM")
# cd("E:\\MyPublications\\MPM_Julia\\Codes\\SpeedTest")

# simple implementation to test the processing speed of data structures versus multi-dimensional arrays

using Printf

mutable struct testDataStructure
    vMember::Array{Float64,2}
end
function testDataStructure()
    testDataStructure(zeros(100,1))
end

function myMain()
iIterations = 10
# Array----------------------------------------------------------------------
# ---------------------------------------------------------------------------
@printf("Creating arrays...")
@time begin
    thisArray01 = Array{Float64, 2}(undef, 1000000,100)
    thisArray02 = Array{Float64, 2}(undef, 1000000,100)
    thisArray03 = Array{Float64, 2}(undef, 1000000,100)
    @printf("...done ")
end
@printf("\n")

@printf("Assigning arrays...")
@time begin
    for j = 1:1:size(thisArray01,2)
        for i = 1:1:size(thisArray01,1)
            thisArray01[i,j] = 1.0*i
            thisArray02[i,j] = -1.0*i
            thisArray03[i,j] = 0
        end
    end
    @printf("...done ")
end
@printf("\n")

@printf("Manipulating arrays ...")
@time begin
    ni = size(thisArray01,1)
    nj = size(thisArray01,2)
    for count in 1:1:iIterations
        for j = 1:1:nj
            for i = 1:1:ni
                thisArray03[i,j] = thisArray01[i,j] + thisArray02[i,j]
            end
        end
    end
    @printf("...done ")
end
@printf("\n")
# println("Array03: ", thisArray03)

# Data Structure-------------------------------------------------------------
# ---------------------------------------------------------------------------
@printf("Creating data structures...")
@time begin
    thisDataStructure01 = Vector{testDataStructure}(undef, 1000000)
    thisDataStructure02 = Vector{testDataStructure}(undef, 1000000)
    thisDataStructure03 = Vector{testDataStructure}(undef, 1000000)
    @printf("...done ")
end
@printf("\n")

@printf("Assigning data structures...")
    @time begin
    for i = 1:1:size(thisArray01,1)
        thisDataStructure01[i] = testDataStructure()
        thisDataStructure02[i] = testDataStructure()
        thisDataStructure03[i] = testDataStructure()
        for j = 1:1:size(thisArray01,2)
            thisDataStructure01[i].vMember[j] = 1.0*i
            thisDataStructure02[i].vMember[j] = -1.0*i
            thisDataStructure03[i].vMember[j] = 0
        end
    end
    @printf("...done ")
end
@printf("\n")

@printf("Manipulating data structure type 1...")
@time begin
    ni = size(thisArray01,1)
    nj = size(thisArray01,2)
    for count in 1:1:iIterations
        for i = 1:1:ni
            for j = 1:1:nj
                thisDataStructure03[i].vMember[j] = thisDataStructure01[i].vMember[j] + thisDataStructure02[i].vMember[j]
            end
        end
    end
    @printf("...done ")
end
@printf("\n")

@printf("Manipulating data structure type 2...")
@time begin
    ni = size(thisArray01,1)
    nj = size(thisArray01,2)
    for count in 1:1:iIterations
        for i = 1:1:ni
            this01 = thisDataStructure01[i]
            this02 = thisDataStructure02[i]
            this03 = thisDataStructure03[i]
            for j = 1:1:nj
                this03.vMember[j] = this01.vMember[j] + this02.vMember[j]
            end
        end
    end
    @printf("...done ")
end
@printf("\n")

@printf("Manipulating data structure type 3...")
@time begin
    ni = size(thisArray01,1)
    nj = size(thisArray01,2)
    for count in 1:1:iIterations
        for i = 1:1:ni
            this01 = thisDataStructure01[i].vMember
            this02 = thisDataStructure02[i].vMember
            this03 = thisDataStructure03[i].vMember
            for j = 1:1:nj
                # thisDataStructure03[i].vMember[j] = thisDataStructure01[i].vMember[j] + thisDataStructure02[i].vMember[j]
                this03[j] = this01[j] + this02[j]
            end
        end
    end
    @printf("...done ")
end
@printf("\n")

# println("Array03: ", thisDataStructure03[2].vMember)
end # myMain

myMain()
