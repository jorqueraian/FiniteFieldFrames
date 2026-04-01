include("Frames.jl")


function maximal_case_O_etf(d::Int, p::Int)
    (d > 1 && gcd(p, d-7) == p) || throw(DomainError((d, p),"p must divide d-7"));

    n = Int(d*(d+1)/2);
    two_sets = collect(combinations(d+1,2));
    gram = zeros(Int, (n, n))

    for verts in combinations(n,2) 
        i = verts[1];
        j = verts[2];
        i_set = two_sets[i];
        j_set = two_sets[j];

        if i_set[1] == j_set[1] || i_set[1] == j_set[2] || i_set[2] == j_set[1] || i_set[2] == j_set[2]
            gram[i,j] = 1;  ## This might be flipped
        else
            gram[i,j] = -1;
        end
    end
    gram += gram'
    matrix(GF(p), gram + 3*I[1:n,1:n])
end

function reconstruct_frame_from_gram(gram, case)
    # see Thm 3.13 and 3.15 of citation needed.
    n = size(gram)[1];
    (n == size(gram)[2]) || throw(DomainError(size(gram),"gram must be square"));

    frame_bool, d = is_frame(gram, case);
    frame_bool || throw(DomainError(size(gram),"gram must be the Gram matrix of a frame"));

    ff = gram.base_ring;

    if case == "O"
        q = 1;
    elseif case == "U"
        (degree(ff) == 2) || throw(DomainError(ff,"in Case U, the provided field must be finite and must be a degree 2 extension"));
        q = Int(sqrt(size(ff)));
    end

    A = matrix(ff, I[1:n,1:n]);
    for j in 2:n 
        for i in 1:(j-1) 
            A[i,j] = gram[i,j] - sum([A[k,i]^q*A[k,j] for k in 1:(i-1)], init=ff(0));
        end
    end
    
    if case == "O"
        for i in 1:n
            num = gram[i,i]- transpose(A[:,i])*A[:,i]
            if !is_square(num)
                Kx, x = ff["x"];
                ff, ffa = finite_field(x^2-num, "a");
                A = ff.(A)
                gram = ff.(gram)
            end
        end
    end
    B = diagonal_matrix(ff(0), n, n);
    
    for i in 1:n 
        if case == "O"
            B[i,i] = sqrt(gram[i,i]- transpose(A[:,i])*A[:,i])
        else
            Kx, x = ff["x"];
            B[i,i] = roots(x^(q+1)-gram[i,i]+conjugate_transpose(A[:,i])*A[:,i])[1]
        end
    end

    Psi = [A; B]
    # Qs: Is Psi always rank = n?, ie is it a basis?
    rank(Psi) == n || throw(DomainError(ff,"frankly i was hopping this wouldnt happen, like every, but it did."));
    if case == "O"
        form = symmetric_form(gram);
        m = symmetric_form(diagonal_matrix([[ff(1) for i in 1:d]; [ff(0) for i in (d+1):n]]));
        #(cong_bool, C) =  is_congruent(m, form);
        # C * gram * C' = daig(1 1 1 1... 0 0 0....)
        # D = C^{-1}
        (cong_bool, D) =  is_congruent(form, m);
        cong_bool || throw(DomainError(ff,"Ok i need to code in the option that you have a symmetric form congruent to the non-square case. I think because we are doing a field extension this wont happen. So if you see this i was wrong."));

        Phi = tranpose(D)[1:d, :]
        (rank(Phi) == d) || throw(DomainError(Phi,"Something broke and I dont know why"));

    elseif case == "U"
        form = hermitian_form(gram);
        m = hermitian_form(diagonal_matrix([[ff(1) for i in 1:d]; [ff(0) for i in (d+1):n]]));
        (cong_bool, D) =  is_congruent(form, m);
        Phi = conjugate_transpose(D)[1:d, :]
        (rank(Phi) == d) || throw(DomainError(Phi,"Something broke and I dont know why"));
    end
    Phi
end
