function B = Bmat(t,T,n)
    tau = pi*t/T;
    fun = 1+(sin(4*tau))/8;
    % fun = 1;
    B = [zeros(n);
         eye(n)*fun];    
end