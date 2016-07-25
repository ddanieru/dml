function L = computeL(M)
    N = (M-1)/2;
    L=zeros(M,M);
    
    for i1=1:M
        i= i1-N-1;
        for j1=1:M
            j=j1-N-1;
            if i==j then
		      L(i1,j1) = 2*%pi*(j^2 + 1);
            end
        end 
    end
    return L;
endfunction
