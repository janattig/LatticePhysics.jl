# COUNTS LINE NUMBERS IN LIBRARY AND UPDATES THEM IN THE MAIN FILE

# list of all line numbers
line_numbers = Int64[]


# open the main file
main_file_filename = "/home/jattig/PackageDevelopement/LatticePhysics.jl/src/LatticePhysics.jl"

# find the folder
main_file_folder = main_file_filename[1:findlast(main_file_filename, '/')]

# get all files within that folder
all_files = readdir(main_file_folder)
# build the absolute filename
#all_files = String["$(main_file_folder)$(f)" for f in all_files]

# open the main file
main_file = open(main_file_filename, "r")
# read the complete main file
main_file_lines = String[rstrip(line) for line in readlines(main_file)]
# close the main file
close(main_file)


# include lines
include_lines = Int64[]

# process all lines
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


# extract the included files
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



# extract line number information
for include_filename in include_files
    # get the number of lines in the corresponding file
    f = open("$(main_file_folder)$(include_filename)", "r")
    push!(line_numbers, length(readlines(f)))
    close(f)
    # print the success
    println("$(include_filename)  ($(line_numbers[end]) lines)")
end

# overwrite the lines in the main file
for i in 1:length(include_lines)
    # get the line
    l = include_lines[i]
    count = line_numbers[i]
    # compose the new line
    line = "# included in subfile ($(count) lines)"
    # substitute the line
    main_file_lines[l-1] = line
end

# overwrite the total count string
main_file_lines[end-2] = "# total lines: $(length(main_file_lines))"
for l in line_numbers
    main_file_lines[end-2] = "$(main_file_lines[end-2]) + $(l)"
end
# add the total number
main_file_lines[end-1] = "# = $(sum(line_numbers) + length(main_file_lines)) lines"
# print the total number
println("  --> total: $(sum(line_numbers) + length(main_file_lines)) lines")

# write to file
open(main_file_filename, "w") do main_file
    for line in main_file_lines
        write(main_file, "$(line)\n")
    end
end
