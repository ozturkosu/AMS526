function [timesL, timesR, fig] = hw1_7(ntimes)
% Template function for testing matrix-matrix and matrix vector
% multiplication.
%
% Place your functions either as separate .m files in the same directory or
% as subfunctions at the end of this .m file.
% INPUTS:
%   NTIMES  - Number of tests to be done. Default value is 5.
% OUTPUTS:
%   TIMESL  - Times (in seconds) for multiplying ABx = (AB)x
%   TIMESR  - Times (in seconds) for multiplying ABx = A(Bx)
%   FIG     - Figure Handle of log-log plot of problem size and total
%             times.
if nargin<1; ntimes = 5; end

% Initialize vectors to store times for multiplication operations as well
% as matrix/vector sizes.
timesL = zeros(ntimes,1);
timesR = timesL;
sizes  = timesR;

% Matrix/vector sizes
sz = 50;
for ii=1:ntimes
    fprintf(1,'sz = %d\n',sz);
    sizes(ii) = sz;
    
    % Initialize matrices Amat and Bmat as sz-by-sz
    Amat = rand(sz);
    Bmat = rand(sz);
    
    % Initialize Xvec as sz-by-1 vector
    Xvec = rand(sz,1);
    
    % Perform ABx = (AB)x
    tic;
    Cmat = MatMatMultiply( Amat, Bmat);
    Yvec = MatVecMultiply( Cmat, Xvec);
    timesL(ii)=toc;
    
    % Perform ABx = A(Bx)
    tic;
    Y1vec = MatVecMultiply( Bmat, Xvec);
    Yvec  = MatVecMultiply( Amat, Y1vec);
    timesR(ii)=toc;
    
    % Double size for next iteration
    sz = sz*2;    
end

% Create Plot
fig = loglog(sizes, timesL, '*-', sizes, timesR, 'r+--');
xlabel('Number of entries');
ylabel('Total Time (sec)');
legend('ABx = (AB)x', 'ABx = A(Bx)', 'Location','NorthWest');

% Print plot into a file hw1_6_result.eps
print -depsc2 hw1_6_result.eps
end

function Yvec = MatVecMultiply(Amat, Xvec)
% Subfunction implementing right matrix-vector multiplication.
assert(size(Amat,2)==size(Xvec,1), 'Matrix and vector sizes are incompatible');
loopSize=size(Xvec,1);
rowNumber=size(Amat,2);
finalSize=size(Xvec,2);
Final=zeros(rowNumber,finalSize);

for i=1:loopSize
    newElement=0;
    for j=1:rowNumber
       newElement=newElement+Amat(i,j)*Xvec(j);
    end
  Final(i)=newElement;
    
end
Yvec=Final;
% TODO: Implement this function

end

function Bmat = MatMatMultiply(A1mat, A2mat)
% Subfunction implementing matrix-matrix multiplication.
assert(size(A1mat,2)==size(A2mat,1), 'Matrices are incompatible for matrix-matrix multiplication.');
rowA1=size(A1mat,1);
colA1=size(A1mat,2);

 
colA2=size(A2mat,2);

FinalMatrix=zeros(rowA1,colA2);
for i1=1:rowA1
     
    for j1=1:colA2
        for k=1:colA1
        FinalMatrix(i1,j1)=FinalMatrix(i1,j1)+(A1mat(i1,k)*A2mat(k,j1));
        end
        
    end
end
Bmat=FinalMatrix;

% TODO: Implement this function
end