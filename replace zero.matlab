function D = convergenceB(M,n)
    d = eps;
    result = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     I = eye(n,n);
     N = zeros(n-1,n-1);
while result ~= 2*n -2
    counter = 0;
%while ~isequal(ismembertol(M,triu(M),eps),ones(size(M)))
    u = M(n,n);
    [Q,R] = qr( M - (u*I) );
    M = R*Q + u*I;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:n-1
                
                
                if abs(M(n,i)) < d  
                     
                    M(n,i)= 0;
                    counter = counter +1;
                end
                if abs(M(i,n)) < d  
                     
                    M(i,n)= 0;
                    counter = counter +1;
                end       

        end
       result = counter;
end
                

       
        D = convergenceA(M,n-1);


    
    
end