# DESIGNER
Diffusion parameter EStImation with Gibbs and NoisE Removal (“DESIGNER”) is an image-processing pipeline for (diffusion) MRI data that has been developed to identify and correct various specific artifacts and confounding factors for an improved accuracy, precision, and robustness in diffusion MRI analysis.  The pipeline optionally corrects for:

1.   Thermal noise
2.   Gibbs ringing
3.   Susceptibility –induced geometric distortions
4.   Eddy current induced spatial distortions and interslice motion
5.   Rician biases using the methods of moments
6.   B1 inhomogeneity
7.   Signal outliers

In general, we recommend using all of the above options, with the exception being that B1 inhomogeneity correction should only be applied if such artifacts are present in the data.

DESIGNER is modular and can easily be tuned to study-specific needs or requirements.  The current implementation is written in Python, but requires the installation of Matlab, MRtrix, and FSL.

The DESIGNER pipeline has been evaluated and described in more technical detail in following publication: Ades-Aron et al. NeuroImage 183: 532-543 (2018).  

Don’t hesitate to reach out to Jelle.Veraart@nyulangone.org or Benjamin.Ades-Aron@nyulangone.org for feedback, suggestions, questions, or comments.
