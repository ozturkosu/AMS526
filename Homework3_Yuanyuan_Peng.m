function Homework3_Yuanyuan_Peng()
    for n = 3:15
    
     m=2*n;
     A=rand(m,n);
 
    
    [Q1,R1]= my_gschmidt(A,n);
    display('Q and R From Gram Schmidt');
    display(Q1);
    display(R1);
   
    [v,R2]=my_householder(A,n);
     display('R and V From Gram Schmidt');
     display(R2);
     display(v);
     %display(B);
    end
   
   
   y= poly(m,n);
   b = impli(v,y,n);
   Q=backsub(y,b);
   display(Q);
end
function y = poly(m,n)
    y1=0;
    for i=1:m 
        ti=(i-1)/(m-1);
    
   
        y1= y1 + i*ti^(i);
        y(i)=y1;
    end
    y= y;
    
end

function b = impli(v,y,n)
m=2*n;
    for k = 1:n
        
       
     
     %    y(k:m)=transpose(y(k:m))-2*v*(transpose(v)*transpose(y(k:m)));
    end
    b=y(k:2*n);
end

function sign_var = my_sign(a)
    if a<0
        sign_var=-1;
    else
        sign_var=1;
    end
end

function [my_Q, my_R] = my_gschmidt(A,n)
    for k = 1:n
        V(:,k)=A(:,k);
    end

   for i = 1:n
       v1=V(:,i);
       R(i,i)=norm(v1);
       Q(:,i)=v1/(R(i,i));
       for j = (i+1):n
           v2=V(:,j);
           q=Q(:,i);
           R(i,j)= transpose(q)*v2;
           V(:,j)=V(:,j)-R(i,j)* Q(:,i);   
       end
   
   end
    my_Q=Q;
    my_R=R;
end
function Q= backsub(y,b)
    Q=y\b;
end
function [v,R] = my_householder(A,n)
        m=2*n;
        for i = 1:n
            X=A(i:m,i);
            x1=X(1);
            
            sign_var=my_sign(x1);
            nor= sign_var*norm(X);
         
            X(1)=X(1)+ nor;
            V(i:m,i)=X;
            V(i:m,i)=V(i:m,i)/norm(V(i:m,i));
            vk=  V(i:m,i);
            A(i:2*n,i:n)=A(i:2*n,i:n)-2*vk*(transpose(vk)*A(i:2*n,i:n)); 
           % b(i:m)=b(i:m)-2*vk*(transpose(vk)*b(i:m));
        end
        R=A;
        v=V;
         
end