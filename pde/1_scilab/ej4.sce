err = zeros(30,1);
for N=1:30

    M = 2*N+1;
    x = zeros(M,1);
    for l=1:M
        x(l)= 2*%pi*l/M;
    end
    u=exp(cos(x));
    V=1+sin(x);
    f=u.*(-sin(x).*sin(x)+cos(x) + V);

    c=solve_nonconstant_coeff(V,f);
    unew=real(fourier_to_real(c));
    err(N)=max(abs(u-unew));
end

plot2d("nl",err);
