# Designer with Matlab Integration

This repo is a collection of scripts aimed at wrapping Designer into MATLAB. This is done for those more familiar with MATLAB than Python. This integration also allows one to use native MATLAB functions such as `for` loops and `dir` to automate processing an entire study.

## What is Designer?

### Original README
The following excerpt is taken from the original designer fork:

>Estimation of diffusion and kurtosis model parameters, including the white matter tract integrity metrics, using Diffusion Kurtosis Imaging
>
>We here provide the code to estimate the diffusion kurtosis tensors from diffusion-weighted images. The (constrained) weighted linear least squares estimator is here preferred because of its accuracy and precision. See “Veraart, J., Sijbers, J., Sunaert, S., Leemans, A. & Jeurissen, B.,  Weighted linear least squares estimation of diffusion MRI parameters: strengths, limitations, and pitfalls. NeuroImage, 2013, 81, 335-346” for more details. Next, a set of diffusion and kurtosis parameter, including the white matter tract integrity metrics, can be calculated form the resulting kurtosis tensor.
>
>Some important notes needs to be considered:
>
>1. Currently, we only provide tools for tensor estimation and parameter calculation. Nonetheless, we recommend that you apply denoising, gibbs correction, motion- and eddy current correction to the diffusion-weighted image prior to tensor fitting. These steps are not included in the provided code, but we are happy to assist (Jelle.Veraart@nyumc.org). 
>
>2. Since the apparent diffusion tensor has 6 independent elements and the kurtosis tensor has 15 elements, there is a total of 21 parameters to be estimated. As an additional degree of freedom is associated with the noise free nondiffusion-weighted signal at least 22 diffusion-weighted images must be acquired for DKI. It can be further shown that there must be at least three distinct b-values, which only differ in the gradient magnitude. Furthermore, at least 15 distinct diffusion (gradient) directions are required (Jensen et al. 2005). Some additional consideration must be made.  The maximal b-value should be chosen carefully and is a trade-off between accuracy and precision. While for DTI, diffusion-weighted images are typically acquired with rather low b-values, about 1000 s⁄mm^2 , somewhat stronger diffusion sensitizing gradients need to be applied for DKI as the quadratic term in the b-value needs to be apparent. It is shown that b-values of about 2000 s⁄mm^2  are sufficient to measure the degree of non-Gaussianity with an acceptable precision (Jensen & Helpern 2010). 
>
>3. Outliers, or “black voxels”, in kurtosis maps are not uncommon. They result from undesired signal fluctuations due to motion, Gibbs ringing, or noise, which can often only be reduced using sophisticated tools.  Unfortunately, those outliers will interfere with the visual and statistical inspection of the kurtosis parameters maps. Smoothing is typically used to suppress those outliers. Use of smoothing must be done with care as image blur partial voluming effects might be introduced. 

### Preprocessing Pipeline
Traditionally, dMRI images are processed with a series of filtering to improve SNR. However, these steps can also muddy the ground truth if one is not careful. The Designer pipleine has been tested to improve SNR wih minimal perturbations to the ground truth.

The table below lists the pipeline differences that between what was traditionally followed at MUSC and Designer

| **MUSC** | **Designer** |
| -------  | ------------ |
| MP-PCA | MP-PCA |
| Gibb's unringing | Gibb's unringing |
| Rician bias correction | Rigid body alignment | 
| Brain mask | EPI + Eddy current correction |
| EPI + Eddy current correction | Brain mask |
| Smooting | Smoothing |
| DKE parameterization | Rician bias correction |
| | DKE parameterization |

The primary advantages Designer has over MUSC's pipeline are rigid body alignment to aid eddy correction, and placement of Rician correction. In our testing of Rician correction placement, it was discovered that placing Rician correction at the end produced significantly different results than when it is placed after Gibb's unringing. According to Dr. Veerart, the Rician bias correction is poorly conditions and depends on the precise and accurate estimation of the Rician distributed signal's expected value. While a greater accuracy in this expectation is achieved perior to motion and Eddy current correction because of signal variance in a local window, additional interpolation steps from said correction smooth the data more to produce a more precise estimation of expected value. This estimation is severely impeded by the presence of noise, as as such, performs better under high SNR. Presently, Dr. Veerart considers the precision gain over potential drops in accuracy due to the interpolation associated to motion and Eddy current distortion correction.

## Installation
This entire software library has the following dependencies:

1. FSL
2. Mrtrix3
3. MATLAB
4. Python

Failing to meet any of these dependencies will prevent Designer from running. As of January 29, 2019, this pipeline was tested and certified to run on MATLAB 2018b and Python 3.6. Due to MATLAB's current support of only up to Python 3.6, version 3.7 has not been tested.

### Install FSL
Obtain and install FSL from their [official page](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/ ). Installation instructions for all compatible operating systems is avaialble on their website.

### Install Mrtrix3
Akin to FSL, installation instructions and downloads are avaialble from Mrtrix3's [official page](http://www.mrtrix.org/). We found homebrew installation to be the most straightforward on a Mac.

### Matlab Installation
MATLAB is very straightforward to install because it's a packaged installer. Just run the installer, wait for it to finish, then proceed to Python installation in the next section.

### Python Installation
There are several python packages available for installation for all operating systems. It is advisable to use a Conda-based distribution such as Miniconda or Anaconda. This guide is based on Miniconda for Ma=c OS Mojave (v.10.14.3), so same directions apart from directory structures will apply.

#### Download Miniconda
Miniconda can be downlaoded from their official page at [Miniconda Latest Releases](https://conda.io/en/latest/miniconda.html). Download and install the most recent release for your system. Version here doesn't matter as we'll be creating a custom Python 3.6 environment within Conda - yes, the flexibility is beautiful.

#### Configure Python 3.6 Environment
Upon successful installation of Miniconda, proceed with the following steps:

Start off by updating conda and then create a new environemnt with a specific version of Python - 3.6 in our case. Name your environment as you please and replace `your_env_name` with your custom name. We used *py36* to keep things simple.

```
conda update conda
conda create -n your_env_name python=3.6
```

Ref: [Managing enviroments (conda)](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

Verify the installation of your environment:
```
conda info --envs
```
Your output should look something like this

![alt-text](https://i.imgur.com/yEqLgoq.png "conda environment locations")

Our environment `py36` appears to have been installed into the environment directory `/miniconda3/envs/py36`. Take note of this directory for the next few steps.

### MATLAB-Python Integration
With the successful set up of MATLAB and Python, we need to ensure that the two languages can communicate with each other. Both languages provide tools to establish some form of communication between each other, and this section will guide you throught he step in ensuring that. It is advisable to perform all configuration via Terminal or command prompt as it minimizes any sources of errors in configuration.

#### Configure Python
Find and located the root of your MATLAB installation by typing the following command in MATLAB's command window:
```
>> matlabroot
ans =

    '/Applications/MATLAB_R2018b.app'
```
Ref: [matlabroot (Mathworks)](https://www.mathworks.com/help/matlab/ref/matlabroot.html)

Take note of this directory, open terminal and activate your custom Python 3.6 environment.
```
source activate your_env_name

Ex:
source activate py36
```
Ref: [Managing environments (conda)](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

In the same terminal window, `cd` to the directory `extern/engines/python` superceding `matlabroot`:
```
cd matlabroot/extern/engines/python

Ex:
cd /Applications/MATLAB_R2018b.app/extern/engines/python
```
This is where you use the installation directory of your custom Python 3.6 environment to run the automatic Python configuration wizard to setup MATLAB Engine API for Python. The `--prefix` flag indicates the directory where `matlab.engine` gets installed.
```
sudo python setup,py install --prefix="your_path_to_conda/envs/your_env_name"

Ex:
sudo python setup.py install --prefix="miniconda3/envs/py36"
```
Ref: [Install MATLAB Engine API for Python in Nondefault Locations](https://www.mathworks.com/help/matlab/matlab_external/install-matlab-engine-api-for-python-in-nondefault-locations.html)

Defauklt installation of Python on Mac creates a symbiolic link to Python dynamic library in `/usr/local/lib`, which does not appear to occur with conda installation. Luckily for us, this can be easily remedied by manually creating a symbiolic link via command line:
```
sudo ln -s path_to_env/lib/libpython/libpython3.6m.dylib /usr/local/lib

Ex:
sudo ln -s /miniconda3/envs/py36/lib/libpython3.6m.dylib /usr/local/lib
```

With these steps, Python should be able to communicate with MATLAB via the matlab.engine package.

#### Configure MATLAB
While it is not necessary to proceed with this step, it is advisable to do this. This ensures that MATLAB's default Python interpretor is is also pointed towards your custom conda environment to minimize any chances of conflict. This is easily achieved by a simple command being directed towards your environment's Python executable.
```
pyversion('/miniconda3/envs/py36/bin/python');
```
Then test whether MATLAB has loaded the libraries correctly by parsing the command pyversion, which should produce an output like below.
```
>> pyversion

       version: '3.6'
    executable: '/miniconda3/envs/py36/bin/python'
       library: '/miniconda3/envs/py36/lib/libpython3.6m.dylib'
          home: '/miniconda3/envs/py36'
      isloaded: 0
```
Once pointed to the correct Python environment, next step is to ensure that MATLAB'S environment path is the same as that defined in system. MATLAB is not aware of any custom paths or environment variables deined in your bash profile. Tho circumvent, this locate your paths and varaibles and add them to MATLAB via the `startup.m` script.

Start off by reading your environment paths into the clipboard via terminal:
```
echo $PATH | pbcopy
```
Then, open MATLAB's `startup.m`. On Mac, this is located in `/Users/your_username/Documents/MATLAB/startup.m`. On any line, add the following command.
```
setenv('PATH','paste/your/clipboard/here')

Ex:
setenv('PATH','/usr/local/Cellar/mrtrix3/3.0_RC3-51-g52a2540d/bin:/Applications/freesurfer/bin:/Applications/freesurfer/fsfast/bin:/Applications/freesurfer/tktools');
```
Ensure that your pasted path is one continous line. Save the file and close it. Restart or start MATLAB, then type the following command to verify your environment paths:
```
getenv('PATH')
```
Do the same for all other variables defined by dependencies (FSL and Mrtrix3) in your bash profile. For Mrtrix3, only a path needs to be defined. FSL, on the other hand, requires a variable FSLDIR that points to it's directory. You can add it to the startup script with:
```
setenv('FSLDIR','path/to/FSL');

Ex:
setenv('FSLDIR','/usr/local/fsl');
```
Once all variables and paths are defined for dependencies, you should have no problems running this software package for automatic study processing. Refer to the screenshots below to to compare our startup.m and ~/.bash_profile.

![alt text](https://imgur.com/efH3T01.png "our ~/.bash_profile")
![a;t text](https://imgur.com/Lkzj1gs.png "out startup.m MATLAB script")

Essentially, try to mimic your startup.m as closely as your bash profile for Mrtrix3 and FSL to ensure a smooth preprocessing pipeline.

Ref: [Calling Python 3 in MATLAB](https://erikreinertsen.com/python3-in-matlab)

### Verify Profile Script
With all binaries and dependencies installed, verify your bash profile script. On Mac, this is `.bash_profile`, which can be accessed by entering the following command in terminal. To make changes, append a sudo at the beginning of this command.
nano ~/.bash_profile
```
nano ~/.bash_profile
```
Ensure that your profile has FSL, MRtrix3 and conda system varibles defined (refer to our .bash_profile's screenshot above) If you do, proceed with a system reboot. Otherwise, check your installation of missing modules.

## Directory Configuration






