function [isolated_real,isolated_nreal,spherical] = q_roots(qpolycoe)
% q_roots computing all zeros of quaternion simple polynomial.
% Implementation of 
% D. Janovsk?a and G. Opfer, ¡°A note on the computation of all zeros of simple quaternionic polynomials,¡±
% SIAM J. Numer. Anal., vol. 48, no. 1, pp. 244¨C256, 2010.

% Input1:  qpolycoe--coefficient vector of quaternion simple polynomial
% (from high order to low order).
% 
% Output1: real isolated zeros.
% Output2: non-real isolated zeros.
% Output3: spherical zeros.


if ~isvector(qpolycoe)
     error('Error. Input must be a vector.')
end

epsl = 10^(-8); % use to chop real numbers that are close to zero by the exact integer 0.


n = length(qpolycoe)-1;  % The order of quaternion polynomial
m = 2*n;  % The order of  its companion polynomial  

%  reverse the coefficient vector such that they are rearranged  from low order to high order.
qpolycoer = reshape(qpolycoe,1,n+1);
qpolycoerev = fliplr(qpolycoer);  
qpolycoecev = reshape(qpolycoerev,n+1,1);


% computing the coefficients of its companion polynomial  
bM1 = repmat(conj(qpolycoerev),n+1,1);
bM2 = repmat(qpolycoecev,1,n+1);
bM = bM1.*bM2;
%bM =scalar(bM);
bM =rot90(bM);

b = zerosq(1,m+1);



for k = 1:m+1
    b(k) = sum(diag(bM,k-n-1));
end

%brev = scalar(fliplr(b));
b = scalar(b);


% computing the zeros of its companion polynomial 
syms x l
 
F(x) = sum(b.*subs(x^l,l,0:m));
S1 = vpasolve(F(x)==0,x);        % vpasolve is more accurate than roots, but is still not accurate enough.

d = double(S1);



d =  reshape(d,1,m);
d(abs(d)<epsl) = 0;


% pick up the real-valued zeros.

zr = d(abs(imag(d))<epsl);
zr = real(zr);
zr = sort(zr);
ind1 = ones(1,length(zr));

for k = 2:length(zr)
    if abs(zr(k)-zr(k-1))<epsl
        ind1(k)=0;
    end
end

zr(ind1==0) = [];

% pick up non-real zeros.

newd = d(imag(d)>epsl);
newd = sort(newd,'ComparisonMethod','abs');
ind2 = ones(1,length(newd));
%newd = unique(newd);

for k = 2:length(newd)
    if abs(newd(k)-newd(k-1))<epsl
        ind2(k)=0;
    end
end

newd(ind2==0) = [];


zi = setdiff(newd, zr);



n1 = length(zi);

zio = zerosq(1,n+1);   % save isolated non-real zeros
zis = zerosq(1,n+1);   % save spherical zeros


% The key step to find isolated non-real spherical zeros, more detail: see
% the reference  D. Janovsk?a and G. Opfer 2010.

for k = 1:n1
    u = real(zi(k))+ 1i*sqrt(abs(zi(k))^2 - real(zi(k))^2);
    alpha = imag(u.^(0:n))./sqrt(abs(zi(k))^2 - real(zi(k))^2);
    beta = zi(k).^(0:n) - alpha*zi(k);
    qalpha = quaternion(real(alpha),imag(alpha),zeros(1,n+1),zeros(1,n+1));
    qbeta = quaternion(real(beta),imag(beta),zeros(1,n+1),zeros(1,n+1));
    A = qalpha*qpolycoecev;
    B = qbeta*qpolycoecev;
    v = conj(A).*B;
    if abs(v)<epsl
        zis(k) =  quaternion(real(zi(k)),imag(zi(k)),0,0);
    else
        x0 = real(zi(k)); 
        x1 = abs(imag(zi(k)))*scalar(v*qi)/abs(vector(v));
        x2 = abs(imag(zi(k)))*scalar(v*qj)/abs(vector(v));
        x3 = abs(imag(zi(k)))*scalar(v*qk)/abs(vector(v));
        x0(abs(x0)<epsl) = 0;
        x1(abs(x1)<epsl) = 0;
        x2(abs(x2)<epsl) = 0;
        x3(abs(x3)<epsl) = 0;
        zio(k) =  quaternion(x0,x1,x2,x3);
    end    
end

zio = zio(abs(zio)>epsl);
zis = zis(abs(zis)>epsl);

% return the different types of zeros.

isolated_real = zr;
isolated_nreal = zio;
spherical = zis;




end

