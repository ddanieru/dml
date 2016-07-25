N = 100;
M = 2*N+1;
x = zeros(M,1);
u = zeros(M,1);
for l=1:M
    x(l)= 2*%pi*l/M;
    //u(l)=1/(2*%pi)*exp(%i*l*x(l));
end
u=1/sqrt(2*%pi)*exp(%i*1*x);
//u=exp(cos(x));
//disp(u);
//plot(x,u)
//V=sin(x);
// Free electron model (Poisson)
//V=x.*0;
// Quantum harmonic oscillator:
k=3
V=0.5*k.*(x-%pi).^2;

//c=solve_nonconstant_coeff(V,f);
[c,lambda]=solve_eig(V);
usol=real(fourier_to_real(c));
disp(inttrap(x,(usol.^2/(2*%pi))),'No. particles:');
//plot(x,usol)
// Wavefunction
plot(x,(usol/sqrt(2*%pi)));
// Probability distribution
plot(x,(usol/sqrt(2*%pi)).^2);



