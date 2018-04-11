function [ diff ] = DiffMat(offDiag,diag,n,m,leftBC,rightBC,dir)
    diff = zeros(n,m);
	
	%{
    for i=1:n
        for j=1:m
            if i==j
                diff(i,j)=diag;
            end
            if (i==j-1 || i==j+1)
                diff(i,j)=offDiag;
            end   
        end
    end
	%}
	
	diff(1,1)=diag;
	if min(n,m)>1
		diff(1,2)=offDiag;
		for i=2:min(n,m)-1
			diff(i,i)=diag;
			diff(i,i+1)=offDiag;
			diff(i,i-1)=offDiag;
		end
		diff(min(n,m),min(n,m))=diag;
		if m<n
			diff(m,m-1)=offDiag;
			diff(m-1,m)=offDiag;
		end
		if n<m
			diff(n,n-1)=offDiag;
			diff(n,n+1)=offDiag;
        end
        if n==m
            diff(m,m-1)=offDiag;
        end
	end
	
    if strcmp(dir,'x')
        j = 2;
        k = 1;
    elseif strcmp(dir,'y')
        j = 1;
        k = 2;
    end
    % BC == 1 => Dirichlet with =0
    if leftBC == 2
        %Reflective
        diff(j,k) = 2*diff(j,k);
    end
    if leftBC == 3
        %Periodic
        if strcmp(dir,'x')
            diff(end,1)= offDiag;
        elseif strcmp(dir,'y')
            diff(1,end)= offDiag;
        end
    end
    
    if strcmp(dir,'x')
        j = n-1;
        k = m;
    elseif strcmp(dir,'y')
        j = n;
        k = m-1;
    end
    % BC == 1 => Dirichlet with =0
    if rightBC == 2
        %Reflective
        diff(j,k) = 2*diff(j,k);
    end
    if rightBC == 3
        %Periodic
        if strcmp(dir,'x')
            diff(1,end)= offDiag;
        elseif strcmp(dir,'y')
            diff(end,1)= offDiag;
        end
    end
end
