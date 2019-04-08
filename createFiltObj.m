function filtObject = createFiltObj(Im, violMask, th, sz)
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
%       Im       --> 3D parameter map or any 3D image
%
%       violMask --> Binary mask containing locations of directional
%                    violations or anywhere a median filter is to be
%                    applied.
%
%       th       --> Percentage threshold for median filtering
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
%       filtObject.OriginalVal
%           - original value of violated voxel in patch
%
%       filtObject.PatchVal
%           - median value of violated voxel in patch
%
%       filtObject.MedianIdx
%           - index of filtObject.PatchVal in a patch sorted in acending
%           order
%
%       filtObject.PatchIdx
%           - index of voxel to use for violation replacement in a patch
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
    % Ensure filter box matrix is odd-sized
elseif mod(sz,2) == 0
    error('Filter matrix size needs to be an odd-numnber');
else
    disp('...applying median filter (1/2)');
end

%% Create Filter
% Distance from centroid to edges of 3D box filter
filtObject.Size = sz;
filtObject.Mask = violMask;
filtObject.PropMask = violMask;
filtObject.Mask(filtObject.Mask >= th) = 1;
filtObject.Mask(filtObject.Mask < th) = 0;
filtObject.Mask = logical(filtObject.Mask);
filtObject.Threshold = th;

centralIdx = median(1:sz);
d2move = abs(sz - centralIdx);

violIdx = find(filtObject.Mask);
filtObject.ViolatedVoxels = numel(violIdx);
filtObject.CorrectedVoxels = 0;
if numel(violIdx) > 0
    disp(sprintf('...%d voxels over threshold being filtered',numel(violIdx)));
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
            filtObject.PatchIdx(i) = NaN;
            continue;
        else
            % Sort all patch values into ascending order and remove NaNs
            patchVals = sort(patchI(find(patchViol == 0)),'ascend');
            patchVals = patchVals(~isnan(patchVals));
            nVals = numel(patchVals);
            
            % Different median algorithms based on whether odd or even
            if mod(nVals,2) == 0    % If even
                % Find the middle two numbers and average them, find the
                % index of patch value closest to mean
                medianIdxTmp = [nVals/2, nVals/2 + 1];
                medMeanTmp = mean([patchVals(medianIdxTmp(1)),...
                    patchVals(medianIdxTmp(2))]);
                
                % Then check if mean is same as patch values from the two
                % central median indexes. If it is, replace it by either
                % medianIdx1 or medianIdx2
                if medMeanTmp == patchVals(medianIdxTmp(1)) & medMeanTmp == patchVals(medianIdxTmp(2))
                    medianIdx = medianIdxTmp(1);
                    
                    % If not, find the index of voxel closest to mean. That
                    % will be the median.
                else
                    medDist = abs([medMeanTmp - patchVals(medianIdxTmp(1)),...
                        medMeanTmp - patchVals(medianIdxTmp(2))]);
                    medDist = single(medDist);
                    % If the distance of mean from two central values is
                    % the same, pick the smaller index medianIdx1
                    if medDist(1) == medDist(2)
                        medianIdx = medianIdxTmp(1);
                    else
                        % Else median index is the one closer to mean
                        medianIdx = medianIdxTmp(find(medDist == min(medDist)));
                    end
                end
                
            else   % If odd
                medianIdx = (nVals + 1) / 2;
            end
            
            filtObject.OriginalVal(i) = Im(I,J,K);
            filtObject.CorrectedVal(i) = patchVals(medianIdx);
            filtObject.MedianIdx(i) = medianIdx;
            % Locate median value in patch and get index. If multiple
            % indexes are present, pick the lowest one
            filtObject.PatchIdx(i) = min(find(...
                patchI == patchVals(medianIdx)));
        end
    catch
        % Some violated voxels may be on the edge of the image and the
        % median filter will fail. This portion will catch such voxels and
        % label the MedianIdx as NaN to specify no median exists. Because
        % the voxel is very close to the edge, a median
        disp(sprintf('Violation voxel [%d,%d,%d] occurs at image edge...skipping',...
            I,J,K));
        filtObject.CorrectedVal(i) = NaN;
        filtObject.MedianIdx(i) = NaN;
        filtObject.PatchIdx(i) = NaN;
    end
    filtObject.X(i) = I;
    filtObject.Y(i) = J;
    filtObject.Z(i) = K;
end
filtObject.CorrectedVoxels = numel(find(~isnan(filtObject.MedianIdx)));
disp(sprintf('...filter object created with %d%% threshold',th*100));
filtObject.FilterStatus = 1;
else
    disp('...No violations found exceeding threshold');
    filtObject.FilterStatus = 0;
end