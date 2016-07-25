function [t1,t2] = test_time_fft(N)
t1=zeros(1,N/5);
t2=zeros(1,N/5);


for i=1:5:N
    [t1(i),t2(i)] = test_fft(N)
end
return t1,t2
endfunction
