function filtObject = createFiltObj(Im, violMask, sz)
%   =======================================================================
%   This function creates a 3D median filter object that specifies the
%   size of the filter, the violation mask, and the index of a median from
%   sorted patch values. The median values and indexes are determined based
%   on the input image, which is usually MK in the case of DKI.
%
%   The output object is required input to the actual median filtering
%   function.
%   -----------------------------------------------------------------------
%   Inputs:
%       I        --> 3D parameter map or any 3D image
%
%       violMask --> Binary mask containing locations of directional
%                    violations or anywhere a median filter is to be
%                    applied.
%
%       sz       --> Size of 3D localized median filter sz x sz x sz
%   -----------------------------------------------------------------------
%   Outputs:
%       filtObject.Size
%           - size of median box filter
%
%       filtObject.Mask
%           - violation mask
%
%       filtObject.ViolatedVoxels
%           - total number of violated voxels defined in violation mask
%
%       filtObject.CorrectedVoxels
%           - total number of violated voxels that were corrected and 'sz'
%           number of pixels away from the image edge
%
%       filtObject.PatchValue
%           - median value of violated voxel patch
%
%       filtObject.MedianIndex
%           - index of filtObject.PatchValue in a patch sorted ina scending
%           order
%   -----------------------------------------------------------------------
%   Author: Siddhartha Dhiman
%   Email:  dhiman@musc.edu
%   Date:   03/22/2019
%   Created with MATLAB 2018b
%   =======================================================================


[Ix Iy Iz] = size(Im);
[Mx My Mz] = size(violMask);
Im = double(Im);

%% Perform Checks
if prod(size(Im)) ~= prod(size(violMask))
    error('Violation mask and parameter map sizes are not equal');
else
    disp('...applying median filter');
end

%% Create Filter
% Distance from centroid to edges of 3D box filter
filtObject.Size = sz;
filtObject.Mask = violMask;

centralIdx = median(1:sz);
d2move = abs(sz - centralIdx);

violIdx = find(violMask);
filtObject.ViolatedVoxels = numel(violIdx);
filtObject.CorrectedVoxels = 0;
for i = 1:length(violIdx);
    
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
            % Sort all patch values into ascending order
            patchVals = sort(patchI(find(patchViol == 0)),'ascend');
            nVals = numel(patchVals);
            
            % Different median algorithms based on whether odd or even
            if mod(nVals,2) == 0    % If even
                % Find the middle two numbers and average them, find the
                % index of patch value closest to mean
                medianIdx = [nVals/2, nVals/2 + 1];
                medMean = mean([patchVals(medianIdx(1)),...
                    patchVals(medianIdx(2))]);
                
                % Then check if mean is same as patch values from the two
                % central median indexes. If it is, replace it by either
                % medianIdx1 or medianIdx2
                if medMean == patchVals(medianIdx(1)) & medMean == patchVals(medianIdx(2))
                    filtObject.PatchValue(i) = patchVals(medianIdx);
                    
                    % If not, find the index of voxel closest to mean. That
                    % will be the median.
                else
                    medDist = abs([medMean - patchVals(medianIdx(1)),...
                        medMean - patchVals(medianIdx(2))]);
                    % If the distance of mean from two central values is
                    % the same, pick the smaller index medianIdx1
                    if medDist(1) == medDist(2)
                        medianIdx = medianIdx(1);
                    else
                        % Else median index is the one closer to mean
                        medianIdx = medianIdx(find(medDist == min(medDist)));
                    end
                    filtObject.PatchValue(i) = patchVals(medianIdx);
                    filtObject.MedianIndex(i) = medianIdx;
                end
            else
                filtObject.PatchValue(i) = patchVals(medianIdx);
                filtObject.MedianIndex(i) = (nVals + 1)/2;
            end
        end
    catch
        disp(sprintf('Violation voxel [%d,%d,%d] occurs at image edge...skipping',...
            I,J,K));
        filtObject.MedianIndex(i) = NaN;
    end
end
filtObject.CorrectedVoxels = numel(find(~isnan(filtObject.MedianIndex)));
end