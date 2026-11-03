function H=r30b_heights(Q)
% Rows: profile cell, then Gauss index. Columns: face cell, then Gauss index.
H=reshape(permute(Q.h,[3,1,4,2]),Q.ns*Q.n,Q.nz*Q.n);
end
