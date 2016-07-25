function c = solve_constant_coeff(f)
    M = size(f,1);
    L = computeL(M);
    fi = 2.*%pi.*real_to_fourier(f);
    c = linsolve(L,-fi);
    //disp(L*c-fi);
    
endfunction
