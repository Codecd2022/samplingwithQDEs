function [coematrix,coePNs] = phicoe(N,s)
% phicoe computes the coefficients of the solution for - ix(k + 1) + jx(k) - ix(k - 1) = x(k)lambda  
% return the coefficient matrix for x(1) to x(N) andthe  coefficient of
% P_N(lambda,s).


% Input-  N: k from 1 to N
%         s: initial value of x(1)

% Output-  coematrix: N by N coefficient matrix for x(1) to x(N)
%          coePNs: coefficient of PNs  % from low to high.

switch nargin
    case 1
       s = 1;
end

CoePhi = zerosq(N+1);
CoePhi(2,2) = s;

for k = 3:N+2
    for n = 2:N+2
        if k>n
            CoePhi(k,n) = quaternion(0,1,0,0)*CoePhi(k-1,n-1) - quaternion(0,0,0,1)*CoePhi(k-1,n) - CoePhi(k-2,n);
        elseif k==n
            CoePhi(k,n) = quaternion(0,1,0,0)*CoePhi(k-1,n-1);
        end
    end
end

coematrix = CoePhi(2:N+1,2:N+1);

coePNs = CoePhi(N+2,2:N+2); % from low to high.

end

