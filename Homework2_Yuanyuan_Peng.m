%Function Condition1
%matrix_gen(3);
for order =2:n
%generate the matrix A
    A=zeros(n);

    for i=1:n
        for j=1:n
         A(i,j)=1/(i+j-1);
        end
    end
    %find A one norm
    
    for k=1:n
        col_sum=1;
        for r=1:n
            col_sum=col_sum+A(r,k);
        
        end
        temp=col_sum;
        
        max=col_sum;
        
    end
    % find L and U and P
    % LU = PA
    [L,U,P]=lu(A);
    
    trans_U = U.'; %U transpose
    trans_L = L.'; %L transpose
    
    %find c with trans_U
    
    c1=ones(n,1); %vecotr with all 1
    c2=c1*-1; %vector with all -1
    %c = -1
    v1=trans_U\c2;
    %c = 1
    v2=trans_U\c1;
    %assigan c value to the c vector and v value to the v vector
    for el=1:n
        if v1>v2
            c(el)=-1;
            v(el)=v1;
        else
            c(el)=1;
            v(el)=v2;
        end
    end
    
    %let L transpose * P equal to some new matrix N
    N=trans_L*P;
    %Ny=v
    y=N\v;
    %Az=y
    z=A\y;
    
    one_z=0;
    for i=1:n
    one_z=one_z+z(i);
    end
    
    one_y=0;
        for j=1:n
        one_y=one_y+y(j);
        end
    A_1_norm1 = one_z/one_y;
    
    
    
    
    
end





%end
 
%function return_value = Condition1(Matr)
  %  disp(Matr);
 %   return_value=1;
%end