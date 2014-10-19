function A = computeSurfAreas(X,tri)
% Computes the area of surface mesh elements.
%
% Input:
% X - vertices
% tri - triangulation
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

A = zeros(size(tri,1),1);
for i=1:length(A)
    v1 = X(tri(i,2),:)-X(tri(i,1),:);
    v2 = X(tri(i,3),:)-X(tri(i,1),:);
    A(i) = 0.5*norm(cross(v1,v2));
end
end