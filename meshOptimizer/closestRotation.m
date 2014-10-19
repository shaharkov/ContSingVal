function [ U,E,V,flip] = closestRotation( A )
% ssvd implementation - return U,E,V so that U*E*V'=A, U,V are rotations, E
% is diagonal with E(end,end)=sign(det(A)), and lastly flip =sign(det(A))
% TODO: this can be improved to use only one determinant computation, of A
% itself as the det of U and V is irrelevant
%
% taken from the code implementing the paper "Injective and Bounded Mappings in 3D".
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%regular svd
[U,E,V]=svd(A);
%compute determinant of V
dV=det(V);
%get source and target dimensions
SD=size(A,2);
TD=size(A,1);
flip=0;
if dV<0 %if V is flipped, simply fix the determinant of both. This still keeps A=UEV'
    V(:,SD)=-V(:,SD);
    U(:,SD)=-U(:,SD);
end
%now det(V)>0, let's see if det(U)<0 - if so det(A)<0 
dU=det(U);
if dU<0
    if TD==SD %if we are of full dimension, only then can there be a flip
        flip=1;
    end
    U(:,TD)=-U(:,TD);%now fix det(U)>0
    E(TD,SD)=-E(TD,SD);%and correct E so that we keep A=UEV'
end
end

