function L = computeLpot(M,V)
    N = (M-1)/2;
    L=zeros(M,M);
    disp(V,"V:");
    for i1=1:M
        i= i1-N-1;
        L(i1,i1) = L(i1,i1) + 2*%pi*(i^2);
        for j1=1:M
            j=j1-N-1;
            if abs(i-j)<N+1 then
               l1=(i-j)+N+1;
               L(i1,j1) = L(i1,j1) + 2*%pi*V(l1);
            end
        end 
    end
    //disp(L);
endfunction


