function u = fourier_to_real(c)
    //M = 2*N+1;
    M= size(c,1);
    N = (M-1)/2;
    u = zeros(M, 1);
    //u = ifft(c);
    for l=1:M
        for k1=1:M
            k = k1 - N - 1;
            u(l) = u(l) + c(k1)*exp(%i*k*2*%pi*l/M);
            //u(l) = u(l)*(-1)^(l+k); 
        end
    end
    //u = fftshift(ifft(ifftshift(c)));
    
    
endfunction
