function sorted_q = q_sort(qvector)
%   sort quaternions in order of  $\succ$ ('descend') defined as follows.
%   For any  $q\in \H$,  $w=w_0+w_1\qi+w_2\qj+w_3\qk,v=v_0+v_1\qi+v_2\qj+v_3\qk\in \theta(q)$,
%   we say $w\succ v$ if there exists $0\leq k\leq 3$ such that $w_k>v_k$ and $w_l=v_l,~ \forall ~0\leq l<k$.


if ~isvector(qvector)
     error('Error. Input must be a vector.')
end

n = length(qvector);

qvector = reshape(qvector,n,1);


X0 = scalar(qvector);
X1 = x(qvector);
X2 = y(qvector);
X3 = z(qvector);

X = [X0,X1,X2,X3];

A = sortrows(X,'descend');

A = transpose(A);



sorted_q  = quaternion(A(1,:),A(2,:),A(3,:),A(4,:));




end

