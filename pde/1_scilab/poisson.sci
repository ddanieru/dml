function [unew,err] = poisson(N,err)

    M = 2*N+1;
    x = zeros(M,1);
    for l=1:M
        x(l)= 2*%pi*(l)/M;
    end
    u=exp(cos(x));
    f=u.*(-sin(x).*sin(x)+cos(x) + 1);

    c=solve_constant_coeff(f);
    disp(c);
    unew=real(fourier_to_real(c));
    err(N)=max(abs(unew-u));
    
endfunction
