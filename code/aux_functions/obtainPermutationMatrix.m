
function Pmn = obtainPermutationMatrix(m,n)

% This function yields a permutation matrix Pmn such that for a matrix A of
% size (mxn) Pmn*vec(A) = vec(A')

E = eye(n);

Pmn = zeros(m*n);
for i = 1:n
    start_id = (i-1)*m+1;
    end_id  = i*m;
    Pmn(:,start_id:end_id) = kron(eye(m),E(:,i));
end

