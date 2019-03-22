function V = localizedMedFilter(I, violMask, sz)
%   This is a function that takes in 3D DTI and DKI maps, along with their
%   directional violation maps and applied a median filter. This function
%   differs from the traditional 3D filter in MATLAB in that it is
%   localized and only applies the filter in location where a violation
%   occurs.
%
%   Inputs:
%       I        --> 3D parameter map or any 3D image
%       violMask --> Binary mask containing locations of directional
%                      violations or anywhere a median filter is to be
%                      applied.
%   Outputs:
%       V --> Output 3D image from the application of median filter.
%
%   Author: Siddhartha Dhiman
%   Date:   03/22/2019
%   Created with MATLAB 2018b


[Ix Iy Iz] = size(I);
[Mx My Mz] = size(violMask);
I = double(I);
V = I;

%% Perform Checks
if any(size(I) ~= size(violMask))
    error('Violation mask and input image size is not the same.')
else
    filterObject = zeros(sz,sz,sz);
end

%% Create Filter
centralIdx = median(1:sz);

% Index of center of square matrix circulating around the perimeter of
% image
iterX = centralIdx:(Ix - (centralIdx-1));
iterY = centralIdx:(Iy - (centralIdx-1));
iterZ = centralIdx:(Iz - (centralIdx-1));

% Nesting 'for' loops here allow a "scanning" like mechanism where the
% region of size sz-by-sz scans from one edge of x-axis to the other, then
% shifts up by one voxel on the y-axis and scans along the x-axis again.
% AFter scanning the entire xy-plane, it moves one voxel up the z-plane and
% rescans the xy-plane. It does this until the entire image is scanned.
for z2 = iterZ
    z1 = z2 - 1;
    z3 = z2 + 1;
    
    for y2 = iterY
        y1 = y2 - 1;
        y3 = y2 + 1;
        
        for x2 = iterX
            x1 = x2 - 1;
            x3 = x2 + 1;
            
            %% Determine Median Values
            
            % Read all map values within this filter
            tmp = I(x1:x3,y1:y3,z1:z3);
            
            % NaN the central value to prevent it from being part of the
            % filter
            tmp(2,2,2) = NaN;
            
            if violMask(x2,y2,z2) == 1
                V(x2,y2,z2) = nanmedian(tmp,'all');
            else
                % do nothing
            end
        end
    end
end
end