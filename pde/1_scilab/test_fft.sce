// test complexity for fourier transform
function [t1,t2] = test_fft(N)


    c=rand(1,N);
    timer();real_to_fourier(fourier_to_real(c));t1=timer();
    timer();fft(ifft(c));t2=timer();
endfunction
