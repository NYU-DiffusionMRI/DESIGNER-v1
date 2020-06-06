classdef rpgdegibbs_s < handle
    
    properties (GetAccess = public, SetAccess = private)
        
    end
    
    methods (Access = public)
        function this = rpgdegibbs_s()
            
        end
        
        function imgdg = degibbs(this,img, dim, pf, params)
            if nargin <= 4
                params = [1 3 20];
            end
            
            [nx,ny,nz,nq] = size(img);
            img = reshape(img,nx,ny,[]);
            
            for k = 1:size(img,3)
                imgi = img(:,:,k);
                if dim==1
                    imgi = imgi.';
                end
                imgi = this.unring(imgi,params);
                if pf==5/8
                    scale = 4;
                    for i = 1:scale
                        imgi(:,i:scale:end) = this.unring_y(imgi(:,i:scale:end),params);
                    end
                elseif pf==6/8
                    scale = 2;
                    for i = 1:scale
                        imgi(:,i:scale:end) = this.unring_y(imgi(:,i:scale:end),params);
                    end
                elseif pf==7/8
                    imgi = imresize(imgi,[nx,ny*3],'nearest');
                    scale = 4;
                    for i = 1:scale
                        imgi(:,i:scale:end) = this.unring_y(imgi(:,i:scale:end),params);
                    end
                    imgi = imresize(imgi,[nx,ny],'nearest');
                elseif pf==1
                else
                    fprintf('Error: rpgdegibbs only supports PF = 5/8, 6/8, and 7/8.');
                end
                
                if dim==1
                    imgi = imgi.';
                end
                img(:,:,k) = imgi;
            end
            imgdg = reshape(img,nx,ny,nz,nq);
        end
        
        function [P, P_pf] = gibbsphantom(this,nP,scale,pf_,sigma)
            nPs = ceil(nP/scale);
            P = phantom(nPs);
            Pk = this.fft2_mri(P);
            scale = (1-scale)/2;
            Pk = Pk(round(end*scale+1):round(end*(1-scale)),round(end*scale+1):round(end*(1-scale)));
            Pk = Pk(1:nP,1:nP);
            Pk = Pk + sigma*randn(nP,nP) + 1j*sigma*randn(nP,nP);
            P = this.ifft2_mri(Pk);
            Pk(:,round(end*pf_+1):end)=0;
            P_pf = this.ifft2_mri(Pk);
        end
        
    end
    
    methods (Static)
        function compilefiles(target)
            mex('-v','-L/usr/local/lib','-lfftw3','-I/usr/local/include/',sprintf('%s/ringRm.cpp',target),'-compatibleArrayDims','-outdir',target);
            mex('-v','-L/usr/local/lib','-lfftw3','-I/usr/local/include/',sprintf('%s/ringRm_y.cpp',target),'-compatibleArrayDims','-outdir',target);
%             mex(sprintf('-v -L/usr/local/lib -lfftw3 -I/usr/local/include/ %s/ringRm_y.cpp -compatibleArrayDims',target));
        end
        
        function v = unring(v,params)
            if nargin == 1
                params = [1 3 20];
            end
            v = double(v);
            v = ringRm(v,params);
        end
        
        function v = unring_y(v,params)
            if nargin == 1
                params = [1 3 20];
            end
            v = double(v);
            v = ringRm_y(v,params);
        end
        
        function imgo = ifft2_mri(imgi)
            n = size(imgi,1)*size(imgi,2);
            imgi = single(imgi); n = single(n);
            imgo = fftshift(fft2(fftshift(imgi)))/sqrt(n);
        end
        
        function imgo = fft2_mri(imgi)
            n=size(imgi,1)*size(imgi,2);
            imgi = single(imgi); n = single(n);
            imgo = fftshift(ifft2(fftshift(imgi)))*sqrt(n);
        end
    end
    
end