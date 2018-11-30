function dwi_obj = unringfun(dwi_obj)
parfor j=1:length(dwi_obj)
    dwi_obj(j).img=unring(dwi_obj(j).img);
    dwi_obj(j).img(isnan(dwi_obj(j).img))=0;
    dwi_obj(j).hdr.dt=[64 0];
end
end