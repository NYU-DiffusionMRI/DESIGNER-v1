function V = applyMedFilt(Im, filtObject)
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
%       Im         --> 3D parameter map or any 3D image
%
%       filtObject --> Median filtering objects created by createFiltObject
%       script. This object specifies values to
%   -----------------------------------------------------------------------
%   Outputs:
%       V          --> output image from median filter
%   -----------------------------------------------------------------------
%   Author: Siddhartha Dhiman
%   Email:  dhiman@musc.edu
%   Date:   03/27/2019
%   Created with MATLAB 2018b
%   =======================================================================

Im = double(Im);

%% Perform Checks
% Check for image and violation map size
if prod(size(Im)) ~= prod(size(filtObject.Mask))
    error('Violation mask and parameter map sizes are not equal');
    % Ensure filter box matrix is odd-sized
elseif mod(filtObject.Size,2) == 0
    error('Filter matrix size needs to be an odd-numnber');
else
    disp('...applying median filter (2/2)');
end

%% Pad Image and Mask
% Apply a nan-padding to all 3 dimensions of image and nan padding to mask; 
% same size as the distance between centroid of patch to edge. This enables
% median filtering of edges.
centralIdx = median(1:filtObject.Size);
d2move = abs(filtObject.Size - centralIdx);
Im = padarray(Im,[d2move d2move d2move],nan,'both');
V = Im;
% violMask = padarray(filtObject.Mask,[d2move d2move d2move],0,'both');

%% Create Filter
if filtObject.FilterStatus == 1
    % Distance from centroid to edges of 3D box filter
%     centralIdx = median(1:filtObject.Size);
%     d2move = abs(filtObject.Size - centralIdx);
%     violIdx = find(violMask);
    
    for i = 1:length(filtObject.PatchIdx)
        
        % If PatchIdx is an integer, perform voxel replacement. Otherwise do
        % nothing and jump to next violation.
        
        if ~isnan(filtObject.PatchIdx(i))
%             [I, J, K] = ind2sub([Ix,Iy,Iz],violIdx(i));
            
            % Index beginning and ending of median filter (box) matrix
            Ib = filtObject.X(i) - d2move;
            Ie = filtObject.X(i) + d2move;
            
            Jb = filtObject.Y(i) - d2move;
            Je = filtObject.Y(i) + d2move;
            
            Kb = filtObject.Z(i) - d2move;
            Ke = filtObject.Z(i) + d2move;
            
%             patchViol = filtObject.Mask(Ib:Ie, Jb:Je, Kb:Ke);
            patchI = Im(Ib:Ie, Jb:Je, Kb:Ke);
            
            V(filtObject.X(i),filtObject.Y(i),filtObject.Z(i)) = ...
                patchI(filtObject.PatchIdx(i));
            
        else
            continue;
        end
    end
else
end
%% Unpad Image
V = V(d2move+1:end-d2move,d2move+1:end-d2move,d2move+1:end-d2move);
end