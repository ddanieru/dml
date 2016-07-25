// Shifting according to:
// http://es.mathworks.com/matlabcentral/newsreader/view_thread/285244

err = zeros(1,1);
for N=1:10
    //N = 10;
    M = 2*N+1;
    x = zeros(M,1);
    for l=1:M
        x(l)= 2*%pi*(l)/M;
    end
    u=exp(cos(x));
    [unew,err] = poisson(N,err);
end


//plot(x,[u,unew])
plot2d("nl",err);


