function coeq = CoeQPoly(q,N,s)
% CoeQPoly Computes the coefficient vector (from low to high) of \langle \vphi(q,s) ,\vphi(\lambda,s)\rangle
%  - ix(k + 1) + jx(k) - ix(k - 1) = x(k)lambda 

% Input-  N: k from 1 to N
%         q: lambda = q
%         s: initial value of x(1)
epsl = 10^(-9); % use to chop real numbers that are close to zero by the exact integer 0.

coematrix = phicoe(N,s);

vphiq = vphi(q,N,s);

MDvPhi = repmat(vphiq,1,N);

Mq =  conj(MDvPhi).*coematrix;

coeq = sum(Mq);
x0 = scalar(coeq);
x1 = x(coeq);
x2 = y(coeq);
x3 = z(coeq);

 x0(abs(x0)<epsl) = 0;
 x1(abs(x1)<epsl) = 0;
 x2(abs(x2)<epsl) = 0;
 x3(abs(x3)<epsl) = 0;

coeq = quaternion(x0,x1,x2,x3);
end

