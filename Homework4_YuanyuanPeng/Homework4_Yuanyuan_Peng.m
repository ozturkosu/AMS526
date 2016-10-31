function Homework4_Yuanyuan_Peng()
% Generate Matrix A and B
A = [2 3 2; 10 3 4; 3 6 1];
B = [6 2 1; 2 3 1; 1 1 1];

display('Testing A');
display(A);
%print out convergence value
r = size(A,1);
matrix_eigA =  QRshift_base(A,r);
display(' QRshift_base');
display(matrix_eigA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Testing B');
display(B);%print out convergence value
c = size(B,1);
matrix_eigB = QRshift(B,c);
display( 'QRshift');
display( matrix_eigB);


end

function C = QRshift_base(M,n) 
  I = eye(n,n);
  M_pre = zeros(n,n);
  while (norm(M-M_pre) >= eps)
    M_pre=M;
    u = M(n,n);
    [Q,R] = qr( M - (u*I) );
    M = R*Q + u*I;
  end
   C = M;
end

function R = QRshift(M,n)
    d = eps;
    result = 0;
    
    I  = eye(n,n);
    N = zeros(n-1,n-1);
    
    while result ~= 2*n-2
        counter = 0;
        u = M(n,n);
        [Q,R] = qr( M - (u*I) );
        M = R*Q + u*I;
        for i = 1:n-1
            if abs(M(n,i)) < d
                M(n,i)=0;
                counter = counter + 1;
            end
            if abs(M(i,n)) < d
                M(i,n)=0;
                counter = counter + 1;
            end
        end
        result =  counter;
    end
    
    for r = 1:n-1
        for c = 1:n-1
            N(r,c)=M(r,c);
        end
    end
    
    Temp = QRshift_base(N,n-1);
    M(1:n-1,1:n-1)=Temp;
    R = M;
end


