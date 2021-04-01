function read_sddp_results(file::AbstractString)
    res_dict = JSON.parsefile(file)
    return convert(Vector{Vector{Dict{String, Any}}}, res_dict)
end
