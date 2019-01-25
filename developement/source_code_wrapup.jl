# COUNTS LINE NUMBERS IN LIBRARY AND UPDATES THEM IN THE MAIN FILE


#
#   HEADER PART, DEFINITION OF MAIN FILE
#

# the main file
main_file_filename = "/home/jattig/PackageDevelopement/LatticePhysics.jl/src/LatticePhysics.jl"

# find the folder
main_file_folder = main_file_filename[1:findlast(main_file_filename, '/')]

# get all files within that folder
all_files = readdir(main_file_folder)

# open the main file
main_file = open(main_file_filename, "r")
# read the complete main file
main_file_lines = String[rstrip(line) for line in readlines(main_file)]
# close the main file
close(main_file)




#
#   EXTRACT NUMBER OF LINES FROM SUBFILES
#

# list of all line numbers
line_numbers = Int64[]

# list of all include lines (line numbers where include(...) is written)
include_lines = Int64[]

# process all lines and find include lines
for l in 5:length(main_file_lines)
    # look at line lines preceeding line l to find include lines
    if  main_file_lines[l-3] == "################################################################################" &&
        main_file_lines[l-2] == "" &&
        length(main_file_lines[l-1]) > 10 &&
        main_file_lines[l-1][1:10] == "# included"
        # this line contains the include
        push!(include_lines, l)
    end
end

# extract the included files from these lines
include_files = String[]
for line in [main_file_lines[i] for i in include_lines]
    # extract the filename
    filename = line[findfirst(line, '"')+1:findlast(line, '"')-1]
    # check if it is in the directory
    if !(filename in all_files)
        println("file not found: $(filename)")
        push!(include_files, "NOTHING.jl")
    else
        # push the filename into the list
        push!(include_files, filename)
    end
end

# extract line number information from the files
for include_filename in include_files
    # get the number of lines in the corresponding file
    f = open("$(main_file_folder)$(include_filename)", "r")
    push!(line_numbers, length(readlines(f)))
    close(f)
    # print the success
    println("$(include_filename)  ($(line_numbers[end]) lines)")
end

# overwrite the lines in the main file specifying the include comments with line number counts
for i in 1:length(include_lines)
    # get the line
    l = include_lines[i]
    count = line_numbers[i]
    # compose the new line
    line = "# included in subfile ($(count) lines)"
    # substitute the line
    main_file_lines[l-1] = line
end

# overwrite the total count string at the end of the file
main_file_lines[end-2] = "# total lines: $(length(main_file_lines))"
for l in line_numbers
    main_file_lines[end-2] = "$(main_file_lines[end-2]) + $(l)"
end

# add the total number
main_file_lines[end-1] = "# = $(sum(line_numbers) + length(main_file_lines)) lines"
# print the total number
println("  --> total: $(sum(line_numbers) + length(main_file_lines)) lines")




#
#   SET COMMENT BLOCKS CORRECTLY
#

# find all comment blocks in main file, included file and specify their data
# main_file_line_start, main_file_line_end, sub_file_lines
comment_blocks = Array{Any,1}[]

# find all comment blocks in subfiles
for (i,include_filename) in enumerate(include_files)
    # get the number of lines in the corresponding file
    f = open("$(main_file_folder)$(include_filename)", "r")
    subfile_lines = readlines(f)
    close(f)
    # find the extent of the comment block
    comment_length = 0
    while length(subfile_lines[comment_length+1]) > 0 && subfile_lines[comment_length+1][1] == '#'
        comment_length += 1
    end
    # add a new entry to the comment block list
    push!(comment_blocks, Any[-1, include_lines[i]-3, subfile_lines[1:comment_length]])
end

# find all comment blocks in main file
for (i,cb) in enumerate(comment_blocks)
    # check the extent
    extent = 0
    while length(main_file_lines[cb[2]-extent]) > 0 && main_file_lines[cb[2]-extent][1] == '#'
        extent += 1
    end
    # save into comment block info
    cb[1] = cb[2] - extent + 1
end

# print all comment blocks
for cb in comment_blocks
    println("from line $(cb[1]) to $(cb[2]):")
    for l in cb[3]
        println(l)
    end
end

# build a new list of all main_file lines
main_file_lines_new = String[]
for l in main_file_lines[1:comment_blocks[1][1]-1]
    push!(main_file_lines_new, l)
end
for i in 1:length(comment_blocks)-1
    for l in comment_blocks[i][3]
        push!(main_file_lines_new, l)
    end
    for l in main_file_lines[comment_blocks[i][2]+1:comment_blocks[i+1][1]-1]
        push!(main_file_lines_new, l)
    end
end
for l in comment_blocks[end][3]
    push!(main_file_lines_new, l)
end
for l in main_file_lines[comment_blocks[end][2]+1:end]
    push!(main_file_lines_new, l)
end

# overwrite original main file lines
main_file_lines = main_file_lines_new

# overwrite the total count string at the end of the file
main_file_lines[end-2] = "# total lines: $(length(main_file_lines))"
for l in line_numbers
    main_file_lines[end-2] = "$(main_file_lines[end-2]) + $(l)"
end

# add the total number
main_file_lines[end-1] = "# = $(sum(line_numbers) + length(main_file_lines)) lines"
# print the total number
println("  --> total: $(sum(line_numbers) + length(main_file_lines)) lines (corrected)")



#
#   LAST STEP: WRITE MAIN FILE AGAIN TO FILE
#

# write to file
open(main_file_filename, "w") do main_file
    for line in main_file_lines
        write(main_file, "$(line)\n")
    end
end
