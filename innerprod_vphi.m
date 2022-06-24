function innerprod = innerprod_vphi(a1,a2,N,s)
%    innerprod_vphi computes inner product of vphi(a1,N,s) and vphi(a2,N,s)
%    related to - ix(k + 1) + jx(k) - ix(k - 1) = x(k)lambda 



innerprod = reshape(conj(vphi(a1,N,s)),1,N)*reshape(vphi(a2,N,s),N,1);





end

