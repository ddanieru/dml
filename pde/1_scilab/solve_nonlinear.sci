function c = solve_nonlinear(V,f,alpha)
    
    eps=10^-6;
    cold = zeros(size(f,1),1);
    for i=1:30
        
        c = solve_nonlin_coeff(V,f);
        c = cold + 0.5*(c-cold);

        unew=fourier_to_real(c);
        V=ones(f)+alpha.*unew.^2;
        fnew=unew.*(-sin(x).*sin(x)+cos(x) + V);
        
        //err = max(abs(f-fnew));
        err = max(abs(c-cold));
        disp(err)
        if err<eps then
            disp(err,'Error=',i,'Number of cycles(converged):')
            break;
        end
        cold=c;
    end
    
    if err>=eps then
        disp(err, 'NOT CONVERGED!! Error is:')
    end


endfunction
