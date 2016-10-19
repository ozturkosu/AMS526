%Function Condition1
%matrix_gen(3);
function Homework2_Yuanyuan_Peng()
tic;
for order =2:12
%order=2;
%generate the matrix A
   A=genMatr(order);
%find A one norm
    %Acond=cond(A,1);
    %disp('cond');
    %disp(Acond);
   A1_norm=Max_col(A);
     
    % find L and U and P
    % LU = PA
    [L,U,P]=lu(A);
    
    trans_U = U.'; %U transpose
    trans_L = L.'; %L transpose
    [V_val, C_val]=findCV(trans_U,order);
    %let some matrix B = trans_L * P use B to find y
     
    Y_Val=findY(trans_L*P,V_val);
    Z_Val=findZ(A, Y_Val);
     
     znorm=norm(Z_Val,1);
     ynorm=norm(Y_Val,1);
    
   
    A_1_norm1 = znorm/ynorm;
    
    mycond = Condition1(A1_norm, A_1_norm1);
     Acond=cond(A,1);
    
   %fig = loglog(order, mycond, 'c*-', order, Acond, 'r+--');
    %plot
    xlabel('order of the matrix');
    ylabel('condition numbers');
    plot(order,mycond,'c*-');
    plot(order,Acond, 'r+--');
    hold on;
    
    %legend('my cond', 'matlab cond', 'Location','NorthWest');
     
end
timesMycond(order)=toc;
disp(timesMycond);
end
%function find z
function zvalue = findZ(Amatr, yVal)
    zvalue=Amatr\yVal;
end
%function find v and c
function [vvalue,cvalue] = findCV(Utrans,orders)
 %find c with trans_U
    n=orders;
    
    c1=ones(n,1); %vecotr with all 1
    c2=c1*-1; %vector with all -1
    v=zeros(n,1);
    %c = -1
    v1=Utrans\c2;
  
    %c = 1
    v2=Utrans\c1;
    
    %assigan c value to the c vector and v value to the v vector
    for el=1:n
        if v1(el)>v2(el)
            c(el)=-1;
            v(el)=v1(el);
        else
            c(el)=1;
            v(el)=v2(el);
        end
    end
    vvalue=v;
    cvalue=c;

end
%function find 
function yvalue =findY(N, Vvalue)
 %Ny=v
    yvalue=N\Vvalue;

end

function Hmatr = genMatr(user_n)
     Hmatr1=zeros(user_n);

    for i=1:user_n
        for j=1:user_n
         Hmatr1(i,j)=1/(i+j-1);
        end
    end
    Hmatr = Hmatr1;
end

function max_col_sum = Max_col(MA)
    n=size(MA,1);
    temp=0;
      max_temp=0;
     for k=1:n
        col_sum=0;
        for r=1:n
            col_sum=col_sum+MA(r,k);
            temp=col_sum;
        end
        if temp>=max_temp
            max_temp=temp;
        end
        
     end
   
    
max_col_sum =max_temp;
end


 
function return_value = Condition1(A1norm, A_1norm)
  
  return_value = A1norm*A_1norm;    
   
end