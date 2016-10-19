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
% function pass in value of A and n 
%TODO finish the funciton with the algo above
	for j = 1:n
		v_row = A(j,:) % assigan jth row of A to the vector v_row 
		for i = 1:(j-1)

		end	
	end
	
end

function result_m2 = my_householder(A,n)
%TODO finish the funciton 
	
	for k = 1:n
		x=Ak:m,k
		vk=sign(x1)norm(x,2)e1+x
		vk = vk/norm(vk,2)
		Ak:m,k:n=Ak:m,k:n-2*vk(vk*Ak:m,k:n)
	end

end
