function [ param_noNaN, param_s ] = repnan( param, mask,akc_out )

  G = fspecial('gaussian',[5 5],0.5);
  param(akc_out==1) = NaN;
  %param(mask==0) = 0;
  xdata = (1:size(param,1))';
  param_i = bsxfun(@(x,y) interp1(y(~isnan(x)),x(~isnan(x)),y),param,xdata);
  %param_i(mask==0) = NaN;
  nslices = size(param,3);
  for s = 1:nslices
      param_s(:,:,s) = nanconv(param_i(:,:,s),G,'nonanout'); 
  end
  %param_s(mask==0) = 0;
  param_noNaN = param; param_noNaN(isnan(param)) = param_s(isnan(param));

end

