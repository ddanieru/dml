N = 40;
M = 2*N+1;
x = zeros(M,1);
for l=1:M
    x(l)= 2*%pi*l/M;
end
u=exp(cos(x));
//disp(u);
//V=1+sin(x);
V = ones(M,1);
alpha=1;
f=u.*(-sin(x).*sin(x)+cos(x) + 1 + alpha.*u.^2);
//f=u.*(-sin(x).*sin(x)+cos(x) + 1);

c= solve_nonlinear(V,f,alpha);
unew=fourier_to_real(c);
//disp(unew)
plot(x,[u,unew])
