from __future__ import division

import os, sys
import numpy as np
import nibabel as nib

def unpatch(X, k):
    kernel = k+k+1
    kx = kernel[0,0]; kx = np.int_(kx)
    ky = kernel[0,1]; ky = np.int_(ky)
    kz = kernel[0,2]; kz = np.int_(kz)
    data = np.zeros((kx, ky, kz, np.int_(X.shape[0])))
    for i in range(0, X.shape[0]):
        data[:,:,:,i] = np.reshape(X[i,:],(kx,ky,kz),order='F')
    return data

def bounding_box(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, cmin, zmin, rmax, cmax, zmax

def sample(mask, kernel):
    k = (kernel-1)/2
    kx = k[0,0]; kx = np.int_(kx)
    ky = k[0,1]; ky = np.int_(ky)
    kz = k[0,2]; kz = np.int_(kz)
    if sampling == 'fast':
        print("Warning: undersampled noise map will be returned")
        bbox = bounding_box(mask)
        n = np.ceil([bbox[3:]/kernel])
        x = np.linspace(np.ceil(bbox[0])+kx+1,np.floor(bbox[0])-kx+bbox[3]+1,n[0,0,0])
        x = np.round(x) - 1
        y = np.linspace(np.ceil(bbox[1])+ky+1,np.floor(bbox[1])-ky+bbox[4]+1,n[0,0,1])
        y = np.round(y) - 1
        z = np.linspace(np.ceil(bbox[2])+kz+1,np.floor(bbox[2])-kz+bbox[5]+1,n[0,0,2])
        z = np.round(z) - 1
        y,x,z = np.meshgrid(x, y, z)
        x = x.flatten(order='F')
        y = y.flatten(order='F')
        z = z.flatten(order='F')
    if sampling == 'full':
        print("Warning: image boundaries are not processed")
        mask[:kx,:,:] = 0
        mask[sx-kx-1:,:,:] = 0
        mask[:,:ky,:] = 0
        mask[:,sy-ky-1:,:] = 0
        mask[:,:,:kz] = 0
        mask[:,:,sz-kz-1:] = 0
           
        x = np.array([], int)
        y = np.array([], int)
        z = np.array([], int)
        for i in range(kz, sz-kz):
           y_,x_ = np.where(mask[:,:,i] == 1)
           x = np.concatenate([x, x_])
           y = np.concatenate([y, y_])
           z = np.concatenate([z, i*np.ones((y_.shape))])
    return x, y, z


################################
################################
# test dataset
datapath = '/Users/Ben/Desktop'
data = os.path.join(datapath,'M0037_034M_DIFF_meso_research.nii')

################################
################################

print('Loading data...')
nii = nib.load(data);
dwi = nii.get_data()

sx, sy, sz, M = dwi.shape
sx = np.int_(sx);
sy = np.int_(sy);
sz = np.int_(sz);
M = np.int_(M);
    
if not 'mask' in globals():
    mask = np.ones((sx,sy,sz))
if not isinstance(mask, bool):
    mask.astype(bool)
    
if not 'kernel' in globals():
    kernel = 5*np.ones((1,3))
if np.isscalar(kernel):
    kernel = kernel*np.ones((1,3))

kernel = kernel+np.mod(kernel,2)-1
k = (kernel-1)/2
kx = k[0,0]; kx = np.int_(kx)
ky = k[0,1]; ky = np.int_(ky)
kz = k[0,2]; kz = np.int_(kz)
N = np.int_(np.prod(kernel))

sampling = 'fast'
if not 'sampling' in globals():
    sampling = 'full'
else:
    sampling = 'fast'

if not 'centering' in globals():
    centering = 0

x, y ,z  = sample(mask, kernel)
xsize = np.int_(x.size)
x = np.int_(x)
y = np.int_(y)
z = np.int_(z)

#Declare variables:
sigma = np.zeros((xsize, 1))
npars = np.zeros((xsize, 1))
signal = np.zeros((M, N, xsize))
    
Sigma = np.zeros((sx, sy, sz))
Npars = np.zeros((sx, sy, sz))
Signal = np.zeros((sx, sy, sz, M))
    
# compute scaling factor for in case N<M
R = np.int_(np.min((M, N)))
scaling = (np.max((M, N))-range(0, R-centering))/N
scaling = scaling[:]

for nn in range(0, 1):
    
    X = dwi[x[nn]-kx:x[nn]+kx+1, y[nn]-ky:y[nn]+ky+1, z[nn]-kz:z[nn]+kz+1, :]
    X = np.reshape(X, (N, M), order='F')
    X = np.transpose(X)
    
    if centering:
        colmean = np.mean(X, axis=1)
        X = X - np.tile(colmean, (M, 1))

    u,vals,v = np.linalg.svd(X, full_matrices=False)
    vals = (vals**2)/N

    csum = np.cumsum(vals[R-centering-1:None:-1])
    cmean = csum[R-centering-1:None:-1]/range(R-centering, 0, -1)
    sigmasq_1 = cmean/scaling
    
    gamma = (M-range(0, R-centering))/N
    rangeMP = 4*np.sqrt(gamma[:])
    rangeData = vals[0:R-centering]-vals[R-centering-1]
    sigmasq_2 = rangeData/rangeMP
    
    t = np.where(sigmasq_2 < sigmasq_1)
    t = t[0][0]

    if not t:
        sigma[nn] = np.nan
        signal[:,:,nn] = X
        t = R+1
    else:
        sigma[nn] = np.sqrt(sigmasq_1[t])
        vals[t:R] = 0
        s = np.matrix(u) * np.diag(np.sqrt(N*vals)) * np.matrix(v)
        if centering:
            s = s + np.tile(colmean, (M, 1))
        
        signal[:,:,nn] = s
    npars[nn] = t

    sys.stdout.write('\r')
    sys.stdout.write("Denoising %f%%" % (100*(nn+1)/xsize))
    sys.stdout.flush()
sys.stdout.write('\n')

for nn in range(0, xsize):
    Sigma[x[nn], y[nn], z[nn]] = sigma[nn]
    Npars[x[nn], y[nn], z[nn]] = npars[nn]
    if sampling == 'fast':
        Signal[x[nn]-kx:x[nn]+kx+1, y[nn]-ky:y[nn]+ky+1, z[nn]-kz:z[nn]+kz+1, :] = unpatch(signal[:,:,nn], k)
    elif sampling == 'full':
        Signal[x[nn], y[nn], z[nn], :] = signal[:, np.int_(np.ceil(N/2))-1, nn]

    sys.stdout.write('\r')
    sys.stdout.write("Unpatching %f%%" % (100*(nn+1)/xsize))
    sys.stdout.flush()
sys.stdout.write('\n')

print('Saving denoised data...')
hdr = nii.header
Signal_ = nib.Nifti1Image(Signal, nii.affine, hdr)
nib.save(Signal_,os.path.join(datapath,'testoutSignal.nii.gz'))
Sigma_ = nib.Nifti1Image(Sigma, nii.affine, hdr)
nib.save(Sigma_,os.path.join(datapath,'testoutSigma.nii.gz'))

