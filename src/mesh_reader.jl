#
# Read mesh data from file
#

using BenchmarkTools



#file = open(mesh_path,"r") 
function getMesh_ASCII(filename)

    mesh_path = mkpath("..//data//"*filename)
    n    = 0
    file    = open("data//mesh.txt","r") 
    n       = 0
    nstart  = 0 
    nend    = 0
    elstart = 0
    eled    = 0

    coord   = Matrix{Float64}[]
    enod    = Array{Int64}[]

    indices = [1 6 7 8 9 10 11]

    # Mark where to reset to
    mark(file)
    for line in readlines(file)
        n += 1
        if chomp(line)=="\$NOD"
            nstart  = n
        elseif chomp(line) == "\$ENDNOD"
            nend    = n
        elseif chomp(line) == "\$ELM"
            elstart = n
        elseif chomp(line) == "\$ENDELM"
            eled    = n
        end
        #println(line)
    end
    reset(file)


    n = 0
    while ! eof(file) 
        n += 1
        line = readline(file)
        if n > nstart + 1 && n < nend - 1
            nums = split(line," ")
            num1 = parse(Float64,String(nums[2]))
            num2 = parse(Float64,String(nums[3]))

            push!(coord, [num1 num2])
        elseif n > elstart + 1 && n < eled - 1
            numb = parse.(Int64, split(line, " "))
            push!(enod,numb[indices])
        else

        end
    end

    close(file)
    return coord, enod
end

#@benchmark getMesh("mesh.txt")
coord, enod = getMesh_ASCII("mesh.txt");


function getEdof(enod)
    edof = Array{Int64,2}(undef,length(enod),2*(size(enod[1],2)-1))
    for el = 1 : enod[end][1]
        edof[el,1:2:11] = enod[el][2:7]*2 - ones(6)
        edof[el,2:2:12] = enod[el][2:7]*2
        println(el)
    end

    return edof
end

edof = getEdof(enod)

function readAscii(filename)
    coord, enod = getMesh_ASCII(filename);
    edof = getEdof(enod)
    return coord, enod, edof
end

filename = "mesh.txt"
 
coord, enod, edof = readAscii(filename)