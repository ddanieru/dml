function c = real_to_fourier(u)
    M = size(u,1);
    N = (M-1)/2;
    c = zeros(M,1);
    //c = fft(u);
    for k1=1:M
        k = k1 - N - 1;
        for l=1:M
            c(k1) = c(k1) + (1/M)*u(l)*exp(-%i*k*2*%pi*l/M);
            //c(k1) = c(k1)*(-1)^(k1+l);
        end
    end
    //c = fftshift(fft(ifftshift(u)));
    
    
endfunction
