function [X, Y, Z] = s2r_ijk2mm(tform, i, j, k)
%S2R_IJK2MM  World coordinates (mm) of voxels given by MATLAB subscripts (1-based).
% niftiinfo(...).Transform expects 0-based indices (NIfTI convention), so 1 is subtracted.
% Check: on the FSL MNI152 2 mm grid, s2r_ijk2mm(tform, 46, 64, 37) returns [0 0 0].
mm = transformPointsForward(tform, [i(:) j(:) k(:)] - 1);
X = mm(:,1);  Y = mm(:,2);  Z = mm(:,3);
end