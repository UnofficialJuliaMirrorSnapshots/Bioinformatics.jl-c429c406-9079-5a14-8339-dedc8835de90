@enum AlignmentType begin
    Global = 1
    Local = 2
    SemiGlobal = 3
end

mutable struct PairwiseAlignment
    seq1::Sequence
    seq2::Sequence
    type::AlignmentType
    sm::Dict{Tuple{Char,Char},Int64}
    d::Int
    align_mat::Array{Float64,2}
    function PairwiseAlignment(seq1, seq2, type, sm, d)
        if type == Global
            align_mat = pairwise_global_alignment_linear_gap(seq1, seq2, sm, d)
            return new(seq1, seq2, type, sm, d, align_mat)
        elseif type == Local
            error("This function is not implemented yet.")
        else
            error("This function is not implemented yet.")
        end
    end
end

"""
    function dotmatrix(s1::Sequence, s2::Sequence)

Calculate the dotplot matrix for given two sequences.
"""
function dotmatrix(s1::Sequence, s2::Sequence)
    mat = zeros(Int8, (length(s1), length(s2)))
    for i in 1:length(s1)
        for j in 1:length(s2)
            if s1[i] == s2[j]
                mat[i, j] = 1
            end
        end
    end
    return mat
end

"""
    function pairwise_global_alignment_linear_gap(seq1, seq2, sm, d)

Needleman-Wunsch algorithm with linear gap penalty.
"""
function pairwise_global_alignment_linear_gap(
    seq1::Sequence,
    seq2::Sequence,
    sm::Dict{Tuple{Char,Char},Int64},
    d::Int
)
    m = length(seq1) + 1
    n = length(seq2) + 1
    mat = zeros(m, n)
    for i in 1:m
        mat[i, 1] = -(i - 1) * d
    end
    for j in 1:n
        mat[1, j] = -(j - 1) * d
    end
    for j in 2:n
        for i in 2:m
            s1 = mat[i-1, j-1] + sm[(seq1[i-1], seq2[j-1])]
            s2 = mat[i-1, j] - d
            s3 = mat[i, j-1] - d
            mat[i, j] = max(s1, s2, s3)
        end
    end
    return mat
end
