function [c,lambda] = solve_eig(V)

    M = size(V,1);
    Vi=real_to_fourier(V);
    L = computeLpot(M,Vi);
    //disp(L)
    //f=u.*(-sin(x).*sin(x)+cos(x) + 1 + sin(x));
    //fi = 2.*%pi.*real_to_fourier(f);
    S=2*%pi*eye(M,M);
    L = (L+L')/2; //trick to make L be symmetric (avoid numerical errors)
//    [lambda,c] = eigs(L/2/%pi,[],1,'SM');
    [V,D] = spec(L/2/%pi); // The solution is the same with L/2*pi and L.
    lambda = D(1,1);
    c = V(:,1);
    disp(lambda,'lowest eigenvalue:')
    //disp(c)
    //disp(L*c)
    //disp(fi)
endfunction
