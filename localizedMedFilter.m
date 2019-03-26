function V = localizedMedFilter(Im, violMask, sz)
%   =======================================================================
%   This is a function that takes in 3D DTI and DKI maps, along with their
%   directional violation maps and applied a median filter. This function
%   differs from the traditional 3D filter in MATLAB in that it is
%   localized and only applies the filter in location where a violation
%   occurs.
%   -----------------------------------------------------------------------
%   Inputs:
%       I        --> 3D parameter map or any 3D image
%       violMask --> Binary mask containing locations of directional
%                    violations or anywhere a median filter is to be
%                    applied.
%       sz       --> Size of 3D localized median filter sz x sz x sz
%   -----------------------------------------------------------------------
%   Outputs:
%       V        --> Output 3D image from the application of median filter.
%   -----------------------------------------------------------------------
%   Author: Siddhartha Dhiman
%   Email:  dhiman@musc.edu
%   Date:   03/22/2019
%   Created with MATLAB 2018b
%   =======================================================================


[Ix Iy Iz] = size(Im);
[Mx My Mz] = size(violMask);
Im = double(Im);
V = Im;

%% Perform Checks
if prod(size(Im)) ~= prod(size(violMask))
    error('Violation mask and parameter map sizes are not equal');
else
    filterObject = zeros(sz,sz,sz);
end

%% Create Filter
% Distance from centroid to edges of 3D box filter
centralIdx = median(1:sz);
d2move = abs(sz - centralIdx);

violIdx = find(violMask);
for i = 1:length(violIdx)
    
    [I, J, K] = ind2sub([Ix,Iy,Iz],violIdx(i));
    
    % Index beginning and ending of median filter (box) matrix
    Ib = I - d2move;
    Ie = I + d2move;
    
    Jb = J - d2move;
    Je = J + d2move;
    
    Kb = K - d2move;
    Ke = K + d2move;
    
    % Place algorithm in a try-catch loop to prevent out-of-index issues
    % from violation pixels too close to edge.
    try
        % Get reference image and violation mask patches
        patchViol = violMask(Ib:Ie, Jb:Je, Kb:Ke);
        patchI = Im(Ib:Ie, Jb:Je, Kb:Ke);
        nViol = numel(find(patchViol));
        
        % Here a check is performed to compute the number of violations in a
        % patch. If all voxels are violations, do nothing. Otherwise, exclude
        % violation voxels from the median calculation
        if nViol == sz^3;
            % If every voxel in patch is a violation, replace nothing
            patchV = patchI;
        else
            % Else exclude all violation voxels from median
            V(I,J,K) = nanmedian(patchI(find(patchViol == 0)),'all');
        end
        
    catch
        disp(sprintf('Violation at [%d,%d,%d] occurs at image edge...skipping',...
            I,J,K));
    end
end
end