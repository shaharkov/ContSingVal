function [ volume ] = computePrimitiveVolume( coords,D )
% compute volume of primitive (tet if D==3, tri if D==2) in the given coords.
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

if D==3
    a=coords(:,1);
    b=coords(:,2);
    c=coords(:,3);
    d=coords(:,4);
    volume=abs(det([a-b b-c c-d]))/6;
else
    a=coords(:,1)-coords(:,3);
    b=coords(:,2)-coords(:,3);
    volume=abs(det([a b]))/2;
end
end

