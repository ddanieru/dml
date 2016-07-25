function c = solve_nonlin_coeff(V,f)
    M = size(f,1);
    Vi=real_to_fourier(V);
    //disp(Vi);
    L = computeLnonl(M,Vi);
    //disp(L)
    fi = 2.*%pi.*real_to_fourier(f);
    c = linsolve(L,-fi);
    //error=L*c-fi;
    finew = L*c;
    fnew = fourier_to_real(finew);
    err = max(abs(real(fnew-f)));

endfunction
