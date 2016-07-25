function c = solve_nonconstant_coeff(V,f)
    M = size(f,1);
    Vi=real_to_fourier(V);
    L = computeLpot(M,Vi);
    fi = 2.*%pi.*real_to_fourier(f);
    c = linsolve(L,-fi);
    //disp(L*c)
    //disp(fi)
endfunction
