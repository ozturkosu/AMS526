for n = 3:15
	A=mx(2*n,n); % A is a matrix

	% call gram schmidt
	%algorithm 7.1 classical gram schmidt
		for j = 1 to n
			vj = aj
			for i = 1 :(j-1)
				rij = qi *qj
				vj = vj - rijqi
			rjj = norm(vj,2)
			qj=vj/rjj
	%
end

function result_m = my_gschmidt(A,n)
% function pass in value of A and n ,
	for j = 1:n
		v_row = A(j,:) % assigan jth row of A to the vector v_row 
		for i = 1:(j-1)
			
		end	
	end
	
end