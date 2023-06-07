function A = Amat(t,T,n)
    tau = pi*t/T;
    fun = (1+cos(2*tau))/4;
    A = [zeros(n) eye(n);
         zeros(n) eye(n)*fun];    
end