using DelimitedFiles

function matrix_from_file(file_loc::String, d::Int, n::Int)
    m = readdlm(file_loc);
    real_mat = reshape(m[1:d*n], (d,n));
    comp_mat = reshape(m[d*n+1:2*d*n], (d,n));
    if iszero(comp_mat)
        real_mat
    else
        real_mat + im*comp_mat
    end
end