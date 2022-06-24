function samples = samplepoint(N,s)
% samplepoint computes sample points associated with
%  - ix(k + 1) + jx(k) -ix(k - 1) = x(k)lambda 


 

 [~,coePNs] = phicoe(N,s); % compute the coefficients of p_N(lambda,s) by phicoe.

 coePNs = fliplr(coePNs);
 
 
 [S1,S2,S3] = q_roots(coePNs); % compute the zeros of p_N(lambda,s) by q_roots.
 
 
 n1 = length(S3);
 n2 = N - length(S1) - length(S2) - length(S3) + 1;
 
 
 newS = zerosq(n1,n2);
 
 % compute sample points by Method 2 described in the paper.
 
for k1 = 1:n1
    q = S3(k1);
    newS(k1,1) = q;
    coeq = CoeQPoly(q,N,s);
    coeq = fliplr(coeq);
    [~,W1,~] = q_roots(coeq);
    W = W1(abs(scalar(W1)-scalar(q)) + abs(abs(vector(W1)) - abs(vector(q)))<10^(-4));
    if isempty(W)
        continue
    end
    W = q_sort(W);
    W2 = W;
    newS(k1,2) = W2(1);
    n3 = length(W2);
    Ind = ones(1,n3);
    for k2 = 2:n3
        T = zerosq(1,k2-1);
        for k3 = 1:k2-1
            T(k3) = innerprod_vphi(W(k3),W2(k2),N,s)*Ind(k3);
        end
        normT = sum(abs(T));
        if normT<10^(-4)
            newS(k1,1+k2) = W2(k2);
        else
            Ind(k2) = 0;
        end
    end
    
end

newS = reshape(newS,1,n1*n2);

newS3 = newS(abs(newS)>10^(-6));


samples = [S1,S2,newS3];


end

