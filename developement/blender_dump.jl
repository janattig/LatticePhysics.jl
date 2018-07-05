# function to dump lattice data to a file for blender to read
function dumpBlenderFile(lattice::Lattice, filename::String; visualize_periodic::Bool=false)
    # open the file in write mode
    f = open(filename, "w")
    # set the radius
    radius = 0.1
    # write all sites
    for p in lattice.positions
        # write start of the line
        write(f, "site:\t")
        # write the position
        if length(p) == 2
            write(f, "$(p[1]), $(p[2]), 0.0")
        elseif length(p) == 3
            write(f, "$(p[1]), $(p[2]), $(p[3])")
        else
            println("site has not correct amount of dimensions!")
        end
        # write the radius
        write(f, ", $(radius)")
        # write line break
        write(f, "\n")
    end
    # close the file
    close(f)
end
