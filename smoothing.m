function dwi = smoothing(dwi, kernelsize, width, csfmask)


% Isotropic smoothing of dwi data, mainly in order to reduce
% the Gibbs ringing. It might be recommended to only smooth the
% high SNR (or low b-valued) data in order not to alter the
% Rice distribution of the low SNR data. This important to maintain
% the high accuracy of the WLLS The max bval can be set (default
% 1000 s/mm^2). The kernelsize is by default [5 5]. The width,
% i.e. the FWHM in voxels, is default 1.25.
% if a csf mask is given as an additional argument, CSF
% infiltration in microstructural signal is avoided during
% smoothing.



    dwi = feval('double', dwi);
    dwi(dwi<=0)=eps;


    if ~exist('width', 'var') || isempty(width)
        width = 1.25;
    end

    if ~exist('kernelsize', 'var') || isempty(kernelsize)
        kernelsize = [5 5];
    elseif numel(kernelsize)==1
        kernelsize = [kernselsize kernelsize];
    end

    if numel(kernelsize)~=2
        kernelsize = [5 5];
    end



    bgmask = isnan(dwi(:,:,:,1));
    for i= 1:size(dwi, 4)

        wmgm = dwi(:,:,:,i); wmgm(csfmask) = NaN;
        wmgm = nansmoothing(wmgm, kernelsize, width);

        csf = dwi(:,:,:,i); csf(~csfmask) = NaN;
        csf = nansmoothing(csf, kernelsize, width);

        total = nansum(cat(4, wmgm, csf), 4);
        total(bgmask) = NaN;

        dwi(:,:,:,i) = total;
    end

end

function B = nansmoothing(A, kernelsize, width)

     if ~exist('kernelsize', 'var') || isempty(kernelsize)
                kernelsize = [5 5];
     elseif numel(kernelsize)==1
                kernelsize = [kernelsize, kernelsize];
     elseif numel(kernelsize)~=2
          warning('incorrect kernel size. Must contain 2 elements. [5 5] selected instead')
          kernelsize = [5 5];
     end
     if ~exist('width', 'var') || isempty(width)
          width = 1.25;
     end


     for i=1:size(A, 3)
        B(:,:,i) = colfilt(A(:,:,i),[5 5],'sliding',@(x)nansmoothing_f(x, kernelsize, width));
     end

     B(isnan(A)) = NaN;
end

function results = nansmoothing_f(segment, kernelsize, width)

    % normal filtering ... NaNs are ignored.
    h = fspecial('gaussian', kernelsize, width/(2*sqrt(2*log(2))));
    H = repmat(h(:), [1 size(segment, 2)]);
    results = nansum(segment.*H);

    % by just ignoring the NaNs the filter no longer adds up to 1! i
    % correct for this by applying a normalization term. i counted the
    % percentage of the filter that has not been applied because of NaNs
    % in the data.

    H_ = H;
    H_(isnan(segment)) = 0;
    n = sum(H_, 1);

    results = results./(n+eps);
end
