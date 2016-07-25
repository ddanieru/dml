//x= [0, %pi/2, %pi, 3*%pi/2, 2*%pi];
M = 5;
x = zeros(M,1);
for l=1:M+1
    x(l)= 2*%pi*(l-1)/M;
end
disp(x)
u=exp(cos(x));
disp(u);
f=u.*(2.*sin(x).*sin(x)-cos(x) + ones(x));
//f=conj(conj(u)'.*(2.*sin(x).*sin(x)-2.*cos(x)))';

c=solve_constant_coeff(f);
disp(fourier_to_real(c));
