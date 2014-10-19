function v = computeVolumes( T,X )
% Computes the volumes/areas for mesh elements.
%
% Input:
% T - triangulation
% X - vertices
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code implementing the paper "Controlling Singular Values with Semidefinite Programming".
% Disclaimer: The code is provided as-is for academic use only and without any guarantees. 
%             Please contact the authors to report any bugs.
% Written by Shahar Kovalsky (http://www.wisdom.weizmann.ac.il/~shaharko/)
%        and Noam Aigerman   (http://www.wisdom.weizmann.ac.il/~noamaig/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(T,2)==3 && size(X,2)==3,
    v = computeSurfAreas(X,T); % surface
else
    v=zeros(size(T,1),1);
    for i=1:size(T,1)
        v(i)=computePrimitiveVolume(X(T(i,:),:)',size(X,2));
    end
end
end

