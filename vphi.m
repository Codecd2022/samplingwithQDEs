function vphiq = vphi(q,N,s)
% vphi computes the column vector (x(1),x(2),...,x(N))^T at lambda=q
%  - ix(k + 1) + jx(k) - ix(k - 1) = x(k)lambda 

% Input-  N: k from 1 to N
%         q: lambda = q
%         s: initial value of x(1)


switch nargin
    case 2
       s = 1;
end

vq = q.^(0:N-1);

vq = reshape(vq,N,1);


coematrix = phicoe(N,s);

vphiq = coematrix*vq;












end

