# Diffusion Kurtosis Imaging
Estimation of diffusion and kurtosis model parameters, including the white matter tract integrity metrics, using Diffusion Kurtosis Imaging

We here provide the code to estimate the diffusion kurtosis tensors from diffusion-weighted images. The (constrained) weighted linear least squares estimator is here preferred because of its accuracy and precision. See “Veraart, J., Sijbers, J., Sunaert, S., Leemans, A. & Jeurissen, B.,  Weighted linear least squares estimation of diffusion MRI parameters: strengths, limitations, and pitfalls. NeuroImage, 2013, 81, 335-346” for more details. Next, a set of diffusion and kurtosis parameter, including the white matter tract integrity metrics, can be calculated form the resulting kurtosis tensor. 

Some important notes needs to be considered:

1. Currently, we only provide tools for tensor estimation and parameter calculation. Nonetheless, we recommend that you apply denoising, gibbs correction, motion- and eddy current correction to the diffusion-weighted image prior to tensor fitting. These steps are not included in the provided code, but we are happy to assist (Jelle.Veraart@nyumc.org). 

2. Since the apparent diffusion tensor has 6 independent elements and the kurtosis tensor has 15 elements, there is a total of 21 parameters to be estimated. As an additional degree of freedom is associated with the noise free nondiffusion-weighted signal at least 22 diffusion-weighted images must be acquired for DKI. It can be further shown that there must be at least three distinct b-values, which only differ in the gradient magnitude. Furthermore, at least 15 distinct diffusion (gradient) directions are required (Jensen et al. 2005). Some additional consideration must be made.  The maximal b-value should be chosen carefully and is a trade-off between accuracy and precision. While for DTI, diffusion-weighted images are typically acquired with rather low b-values, about 1000 s⁄mm^2 , somewhat stronger diffusion sensitizing gradients need to be applied for DKI as the quadratic term in the b-value needs to be apparent. It is shown that b-values of about 2000 s⁄mm^2  are sufficient to measure the degree of non-Gaussianity with an acceptable precision (Jensen & Helpern 2010). 

3. Outliers, or “black voxels”, in kurtosis maps are not uncommon. They result from undesired signal fluctuations due to motion, Gibbs ringing, or noise, which can often only be reduced using sophisticated tools.  Unfortunately, those outliers will interfere with the visual and statistical inspection of the kurtosis parameters maps. Smoothing is typically used to suppress those outliers. Use of smoothing must be done with care as image blur partial voluming effects might be introduced. 


DESIGNER.py is a python script that performs complete diffusion weighted image processing. It uses libraries from FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki), mrtrix3 (http://www.mrtrix.org/), and MATLAB to perform DWI denoising, Gibbs artifact correction, Rician bias correction, EPI distortion correction, Eddy current correction, motion correction, B1 bias field correction, outlier rejection, and Tensor fitting. This script will run on either Linux or Mac OS X as long as the shell environement is set up correctly. 
Users interested in running DESIGNER should follow the following steps:
1.  download DESIGNER.py as well as the supporting MATLAB functions above. We also rely on Matlab tools for Nifti and Analyze Images (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
2.  Set up the python engine for matlab: 

    cd "matlabroot/extern/engines/python"

    python setup.py install

3.  Make sure that mrtrix3 and fsl binaries are correctly sourced in your shell environment PATH viarable:

    MRTRIXDIR=<path/to/mrtrix3>

    FSLDIR=<path/to/fsl>

    . ${FSLDIR}/etc/fslconf/fsl.sh

    PATH=${FSLDIR}/bin:${PATH}

    PATH=${MRTRIXDIR}/bin:${PATH}

    export FSLDIR PATH

4.  Gibbs correction is currently implemented using "unring" (https://bitbucket.org/reisert/unring). This function must be downloaded and compiled using either the matlab mex-file or using fsl-source libraries, the location of the unring executable you want to use (either matlab or fsl) should also be added to the environment PATH variable:

    UNRING=<path/to/unring>

    depending on whether you compile unring using matlab or fsl, do:

    PATH=${UNRING}/matlab:${PATH}

    or

    PATH=${UNRING}/fsl:${PATH}

