#!/usr/bin/env python

# script that runs DESIGNER
import matlab.engine
import os
PATH = os.environ['PATH'].split(":")
mrtrixbin = [s for s in PATH if "mrtrix3" in s][0]
if not mrtrixbin:
    print("cannot find path to mrtrix3, please make sure <path/to/mrtrix3/bin> is in your PATH")
    quit()
mrtrixlib = "".join(mrtrixbin)[:-3]+'lib'
print(mrtrixlib)

import inspect, sys, numpy as np, math, gzip, shutil
from distutils.spawn import find_executable
sys.path.insert(0, mrtrixlib)
from mrtrix3 import app, file, fsl, image, path, phaseEncoding, run

app.init('Benjamin Ades-Aron (benjamin.ades-aron@nyumc.org)',
	'DWI processing with DESIGNER')
app.cmdline.addDescription("""1. pre-check: concatenate all dwi series and make sure the all diffusion AP images and PA image have the same matrix dimensions/orientations 
 						2. denoise the complete dataset\n
                        
						3. gibbs ringing correction on complete dataset\n
                        
 						4. If multiple diffusion series are input, do rigid boy alignment\n

 						5. topup + eddy, rotate output bvecs\n
                        
 						6. perform b1 bias correction on each dwi subset\n
                        
 						7. CSF excluded smoothing\n
                        
                        8. Rician bias correction\n
                        
                        9. Normalization to white matter in first b0 image\n
                        
 						10. irwlls outlier map, cwlls dki fit\n
                        
 						11. outlier detection and removal
 	""")
app.cmdline.addCitation('','Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-aron, B.; Sijbers, J. & Fieremans, E. Denoising of diffusion MRI using random matrix theory. NeuroImage, 2016, 142, 394-406, doi: 10.1016/j.neuroimage.2016.08.016',True)
app.cmdline.addCitation('','Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping using random matrix theory. Magn. Res. Med., 2016, 76(5), 1582-1593, doi:10.1002/mrm.26059',True)
app.cmdline.addCitation('','Kellner, E., et al., Gibbs-Ringing Artifact Removal Based on Local Subvoxel-Shifts. Magnetic Resonance in Medicine, 2016. 76(5): p. 1574-1581.',True)
app.cmdline.addCitation('','Koay, C.G. and P.J. Basser, Analytically exact correction scheme for signal extraction from noisy magnitude MR signals. Journal of Magnetic Resonance, 2006. 179(2): p. 317-322.',True)
app.cmdline.addCitation('', 'Andersson, J. L. & Sotiropoulos, S. N. An integrated approach to correction for off-resonance effects and subject movement in diffusion MR imaging. NeuroImage, 2015, 125, 1063-1078', True)
app.cmdline.addCitation('', 'Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister, P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; De Stefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004, 23, S208-S219', True)
app.cmdline.addCitation('', 'Skare, S. & Bammer, R. Jacobian weighting of distortion corrected EPI data. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 5063', True)
app.cmdline.addCitation('', 'Andersson, J. L.; Skare, S. & Ashburner, J. How to correct susceptibility distortions in spin-echo echo-planar images: application to diffusion tensor imaging. NeuroImage, 2003, 20, 870-888', True)
app.cmdline.addCitation('','Zhang, Y.; Brady, M. & Smith, S. Segmentation of brain MR images through a hidden Markov random field model and the expectation-maximization algorithm. IEEE Transactions on Medical Imaging, 2001, 20, 45-57',True)
app.cmdline.addCitation('', 'Smith, S. M.; Jenkinson, M.; Woolrich, M. W.; Beckmann, C. F.; Behrens, T. E.; Johansen-Berg, H.; Bannister P. R.; De Luca, M.; Drobnjak, I.; Flitney, D. E.; Niazy, R. K.; Saunders, J.; Vickers, J.; Zhang, Y.; DeStefano, N.; Brady, J. M. & Matthews, P. M. Advances in functional and structural MR image analysis and implementation as FSL. NeuroImage, 2004,23, S208-S219',True)
app.cmdline.addCitation('','Collier, Q., et al., Iterative reweighted linear least squares for accurate, fast, and robust estimation of diffusion magnetic resonance parameters. Magn Reson Med, 2015. 73(6): p. 2174-84.',True)
app.cmdline.add_argument('input',  help='The input DWI series. For multiple input series, separate file names with commas (i.e. dwi1.nii,dwi2.nii,...)')
app.cmdline.add_argument('output', help='The output directory (includes diffusion parameters, kurtosis parameters and processed dwi) unless option -processing_only is used, in which case this is the output basename')
options = app.cmdline.add_argument_group('Other options for the DESIGNER script')
options.add_argument('-denoise', action='store_true', help='Perform dwidenoise')
options.add_argument('-extent', metavar=('<size>'), help='Denoising extent. Default is 5,5,5')
options.add_argument('-degibbs', action='store_true', help='Perform Gibbs artifact correction')
options.add_argument('-rician', action='store_true', help='Perform Rician bias correction')
options.add_argument('-rician_lowsnr', metavar=('<b-value where SNR=2>'), help='Perform Rician bias correction specifying the b-value where snr<2')
options.add_argument('-prealign', action='store_true', help='If there are multiple input diffusion series, do rigid alignment prior to eddy to maximize motion correction performance')
options.add_argument('-eddy', action='store_true', help='run fsl eddy (note that if you choose this command you must also choose a phase encoding option')
options.add_argument('-b1correct', action='store_true', help='Include a bias correction step in dwi preprocessing', default=False)
options.add_argument('-normalise', action='store_true', help='normalize the dwi volume to median b0 CSF intensity of 1000 (useful for multiple dwi acquisitions)', default=False)
options.add_argument('-smooth', metavar=('<fwhm>'), help='Include a (csf-free) smoothing step during dwi preprocessing. FWHM is usually 1.20 x voxel size')
options.add_argument('-DTIparams', action='store_true', help='Include DTI parameters in output folder (md,ad,rd,fa,eigenvalues, eigenvectors')
options.add_argument('-DKIparams', action='store_true', help='Include DKI parameters in output folder (mk,ak,rk)')
options.add_argument('-WMTIparams', action='store_true', help='Include DKI parameters in output folder (awf,ias_params,eas_params)')
options.add_argument('-akc', action='store_true', help='brute force K tensor outlier rejection')
options.add_argument('-kcumulants', action='store_true',help='output the kurtosis tensor with W cumulant rather than K')
options.add_argument('-mask', action='store_true',help='compute a brain mask prior to tensor fitting to stip skull and improve efficientcy')
options.add_argument('-datatype', metavar=('<spec>'), help='If using the "-processing_only" option, you can specify the output datatype. Valid options are float32, float32le, float32be, float64, float64le, float64be, int64, uint64, int64le, uint64le, int64be, uint64be, int32, uint32, int32le, uint32le, int32be, uint32be, int16, uint16, int16le, uint16le, int16be, uint16be, cfloat32, cfloat32le, cfloat32be, cfloat64, cfloat64le, cfloat64be, int8, uint8, bit')
options.add_argument('-fit_constraints',help='constrain the wlls fit (default 0,1,0)')
options.add_argument('-outliers',action='store_true',help='Perform IRWLLS outlier detection')
options.add_argument('-fslbvec',metavar=('<bvecs>'),help='specify bvec path if path is different from the path to the dwi or the file has an unusual extention')
options.add_argument('-fslbval',metavar=('<bvals>'),help='specify bval path if path is different from the path to the dwi or the file has an unusual extention')
rpe_options = app.cmdline.add_argument_group('Options for specifying the acquisition phase-encoding design')
rpe_options.add_argument('-rpe_none', action='store_true', help='Specify that no reversed phase-encoding image data is being provided; eddy will perform eddy current and motion correction only')
rpe_options.add_argument('-rpe_pair', metavar=('<reverse PE b=0 image>'), help='Specify the reverse phase encoding image')
rpe_options.add_argument('-rpe_all', metavar=('<reverse PE dwi volume>'), help='Specify that ALL DWIs have been acquired with opposing phase-encoding; this information will be used to perform a recombination of image volumes (each pair of volumes with the same b-vector but different phase encoding directions will be combined together into a single volume). The argument to this option is the set of volumes with reverse phase encoding but the same b-vectors as the input image')
rpe_options.add_argument('-rpe_header', action='store_true', help='Specify that the phase-encoding information can be found in the image header(s), and that this is the information that the script should use')
rpe_options.add_argument('-pe_dir', metavar=('<phase encoding direction>'), help='Specify the phase encoding direction of the input series (required if using the eddy option). Can be a signed axis number (e.g. -0, 1, +2), an axis designator (e.g. RL, PA, IS), or NIfTI axis codes (e.g. i-, j, k)')
app.parse()

def splitext_(path):
    for ext in ['.tar.gz', '.tar.bz2','.nii.gz']:
        if path.endswith(ext):
            return path[:-len(ext)], path[-len(ext):]
    return os.path.splitext(path)

designer_root = os.path.dirname(os.path.realpath(__file__))
DKI_root = os.path.abspath(os.path.join(designer_root,'..'))

app.makeTempDir()

fsl_suffix = fsl.suffix()

UserCpath = app.args.input.rsplit(',')
DWIlist = [os.path.realpath(i) for i in UserCpath]

isdicom = False
for i in DWIlist:
    if not os.path.exists(i):
        print('cannot find input ' + i)
        quit()
    if os.path.isdir(i):
        format=image.Header(i).format()
        if format == 'DICOM':
            isdicom = True
        else:
            print('input is a directory but does not contain DICOMs, quitting')
            quit()

DWIflist = [splitext_(i) for i in DWIlist]
DWInlist = [i[0] for i in DWIflist]
DWIext = [i[1] for i in DWIflist]
miflist = []
idxlist = []
dwi_ind_size = [[0,0,0,0]]

if not app.args.fslbval:
    bvallist = [i + '.bval' for i in DWInlist]
else:
    UserBvalpath = app.args.fslbval.rsplit(',')
    bvallist = [os.path.realpath(i) for i in UserBvalpath]
if not app.args.fslbvec:
    bveclist = [i + '.bvec' for i in DWInlist]
else:
    UserBvecpath = app.args.fslbvec.rsplit(',')
    bveclist = [os.path.realpath(i) for i in UserBvecpath]

if len(DWInlist) == 1:
    if not isdicom:
        run.command('mrconvert -stride -1,2,3,4 -fslgrad ' + bveclist[0] + ' ' + bvallist[0] + ' ' + ''.join(DWInlist) + ''.join(DWIext) + ' ' + path.toTemp('dwi.mif',True))
    else:
        run.command('mrconvert -stride -1,2,3,4 ' + ''.join(DWInlist) + ' ' + path.toTemp('dwi.mif',True))
else:
    for idx,i in enumerate(DWInlist):
        if not isdicom:
            run.command('mrconvert -stride -1,2,3,4 -fslgrad ' + bveclist[idx] + ' ' + bvallist[idx] + ' ' + i + DWIext[idx] + ' ' + path.toTemp('dwi' + str(idx) + '.mif',True))
        else:
            run.command('mrconvert -stride -1,2,3,4 ' + i + ' ' + path.toTemp('dwi' + str(idx) + '.mif',True))
        dwi_header = image.Header(path.toTemp('dwi' + str(idx) + '.mif',True))
        dwi_ind_size.append([ int(s) for s in dwi_header.size() ])
        miflist.append(path.toTemp('dwi' + str(idx) + '.mif',True))
    DWImif = ' '.join(miflist)
    run.command('mrcat -axis 3 ' + DWImif + ' ' + path.toTemp('dwi.mif',True))

app.gotoTempDir()

# get diffusion header info - check to make sure all values are valid for processing
dwi_header = image.Header(path.toTemp('dwi.mif',True))
dwi_size = [ int(s) for s in dwi_header.size() ]
grad = dwi_header.keyval()['dw_scheme']
grad = [ line for line in grad ]
grad = [ [ float(f) for f in line ] for line in grad ]
stride = dwi_header.strides()
num_volumes = 1
if len(dwi_size) == 4:
  num_volumes = dwi_size[3]
bval = [i[3] for i in grad]

nvols = [i[3] for i in dwi_ind_size]
for idx,i in enumerate(DWInlist):
    if len(DWInlist) == 1:
        tmpidxlist = range(0,num_volumes)
    else:
        tmpidxlist = range(sum(nvols[:idx+1]),sum(nvols[:idx+1])+nvols[idx+1])
    idxlist.append(','.join(str(i) for i in tmpidxlist))

# Perform initial checks on input images
if not grad:
  app.error('No diffusion gradient table found')
if not len(grad) == num_volumes:
  app.error('Number of lines in gradient table (' + str(len(grad)) + ') does not match input image (' + str(num_volumes) + ' volumes); check your input data')

if app.args.extent:
	extent = app.args.extent
else: extent = '5,5,5'

run.command('mrconvert dwi.mif working.mif')

# denoising
if app.args.denoise:
    print("...Beginning denoising")
    run.command('dwidenoise -extent ' + extent + ' -noise fullnoisemap.mif working.mif dwidn.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwidn.mif working.mif')

# gibbs artifact correction
if app.args.degibbs:
    print("...Beginning degibbsing")
    run.command('mrdegibbs -nshifts 20 -minW 1 -maxW 3 working.mif dwigc.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwigc.mif working.mif')

# pre-eddy alignment for multiple input series
if app.args.prealign:
    if len(DWInlist) != 1:
        miflist = []
        for idx,i in enumerate(DWInlist):
            run.command('mrconvert -coord 3 ' + idxlist[idx] + ' working.mif dwipretf' + str(idx) + '.mif')
            run.command('dwiextract -bzero dwipretf' + str(idx) + '.mif - | mrconvert -coord 3 0 - b0pretf' + str(idx) + '.mif')
            if idx > 0:
                run.command('mrregister -type rigid -noreorientation -rigid rigidXform' + str(idx) + 'to0.txt b0pretf' + str(idx) + '.mif b0tf0.mif')
                run.command('mrtransform -linear rigidXform' + str(idx) + 'to0.txt dwipretf' + str(idx) + '.mif dwitf' + str(idx) + '.mif')
                miflist.append('dwitf' + str(idx) + '.mif')
        DWImif = ' '.join(miflist)
        run.command('mrcat -axis 3 dwitf0.mif ' + DWImif + ' dwitf.mif')
        run.function(os.remove,'working.mif')
        run.command('mrconvert dwitf.mif working.mif')

# epi + eddy current and motion correction
# if number of input volumes is greater than 1, make a new acqp and index file.
if app.args.eddy:
    print("...Beginning EDDY")
    if app.args.rpe_none:
        run.command('dwipreproc -eddy_options " --repol --data_is_shelled" -rpe_none -pe_dir ' + app.args.pe_dir + ' working.mif dwiec.mif')
    elif app.args.rpe_pair:
        run.command('dwiextract -bzero dwi.mif - | mrconvert -coord 3 0 - b0pe.mif')
        rpe_size = [ int(s) for s in image.Header(path.fromUser(app.args.rpe_pair,True)).size() ]
        if len(rpe_size) == 4:
            run.command('mrconvert -coord 3 0 ' + path.fromUser(app.args.rpe_pair,True) + ' b0rpe.mif')
        else: run.command('mrconvert ' + path.fromUser(app.args.rpe_pair,True) + ' b0rpe.mif')
        run.command('mrcat -axis 3 b0pe.mif b0rpe.mif rpepair.mif')
        run.command('dwipreproc -eddy_options " --repol --data_is_shelled" -rpe_pair -se_epi rpepair.mif -pe_dir ' + app.args.pe_dir + ' working.mif dwiec.mif')
    elif app.args.rpe_all:
        run.command('mrconvert -export_grad_mrtrix grad.txt dwi.mif tmp.mif')
        run.command('mrconvert -grad grad.txt ' + path.fromUser(app.args.rpe_all,True) + ' dwirpe.mif')
        run.command('mrcat -axis 3 working.mif dwirpe.mif dwipe_rpe.mif')
        run.command('dwipreproc -eddy_options " --repol --data_is_shelled" -rpe_all -pe_dir ' + app.args.pe_dir + ' dwipe_rpe.mif dwiec.mif')
        run.function(os.remove,'tmp.mif')
    elif app.args.rpe_header:
        run.command('dwipreproc -eddy_options " --repol --data_is_shelled" -rpe_header dwipe_rpe.mif dwiec.mif')
    elif not app.args.rpe_header and not app.args.rpe_all and not app.args.rpe_pair:
        print("the eddy option must run alongside -rpe_header, -rpe_all, or -rpe_pair option")
        quit()
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwiec.mif working.mif')

# b1 bias field correction
if app.args.b1correct:
    print("...Beginning B1 correction")
    if len(DWInlist) == 1:
        run.command('dwibiascorrect -fsl working.mif dwibc.mif')
    else:
        # b1 correction may still need to be done individually for each diffusion series ...
        miflist = []
        for idx,i in enumerate(DWInlist):
            run.command('mrconvert -coord 3 ' + idxlist[idx] + ' working.mif dwiprebc' + str(idx) + '.mif')
            run.command('dwibiascorrect -fsl dwiprebc' + str(idx) + '.mif dwibc' + str(idx) + '.mif')
            miflist.append('dwibc' + str(idx) + '.mif')
            DWImif = ' '.join(miflist)
        run.command('mrcat -axis 3 ' + DWImif + ' dwibc.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwibc.mif working.mif')

# generate a final brainmask
if app.args.mask or app.args.smooth or app.args.normalise:
    print("...Computing brain mask")
    run.command('dwiextract -bzero working.mif - | mrmath -axis 3 - mean b0bc.nii')
    # run.command('dwi2mask dwibc.mif - | maskfilter - dilate brain_mask.nii')
    # run.command('fslmaths b0bc.nii -mas brain_mask.nii brain')
    run.command('bet b0bc.nii brain' + fsl_suffix + ' -m -f 0.25')
if os.path.isfile('brain_mask.nii.gz'):
    with gzip.open('brain_mask' + fsl_suffix, 'rb') as f_in, open('brain_mask.nii', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    run.function(os.remove,'brain_mask' + fsl_suffix)

if app.args.smooth or app.args.normalise:
    print("...Computing CSF mask")
    run.command('mrconvert -force -export_grad_mrtrix grad.txt working.mif dwibc.nii')
    run.command('fast -n 4 -t 2 -o tissue brain' + fsl_suffix)
    csfclass = []
    for i in range(4):
        run.command('fslmaths tissue_pve_' + str(i) + fsl_suffix + ' -thr 0.95 -bin tissue_pve_thr' + str(i) + fsl_suffix)
        csfclass.append(float(run.command('fslstats brain' + fsl_suffix + ' -k ' + 'tissue_pve_thr' + str(i) + fsl_suffix + ' -P 95')[0]))
    csfind = np.argmax(csfclass)
    run.command('fslmaths tissue_pve_' + str(csfind) + fsl_suffix + ' -thr 0.7 -bin CSFmask' + fsl_suffix)
    if fsl_suffix.endswith('.nii.gz'):
        run.command('mrconvert CSFmask' + fsl_suffix + ' CSFmask.nii')

# smoothing (excluding CSF)
if app.args.smooth:
    print("...Beginning smoothing")
    os.chdir(designer_root)
    eng = matlab.engine.start_matlab()
    eng.runsmoothing(path.toTemp('',True),app.args.smooth,DKI_root,nargout=0)
    eng.quit()
    app.gotoTempDir()
    run.command('mrconvert -grad grad.txt dwism.nii dwism.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwism.mif working.mif')

run.command('mrinfo -export_grad_fsl dwi_designer.bvec dwi_designer.bval working.mif')

# rician bias correction
if app.args.rician and app.args.rician_lowsnr:
    print("...Choosing low snr Rician Bias correction")
    app.args.rician = 0
if app.args.rician:
    print("...Beginning Rician correction")
    if app.args.denoise:
        run.command('mrcalc fullnoisemap.mif -finite fullnoisemap.mif 0 -if lowbnoisemap.mif')
        run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
    else:
        run.command('dwidenoise -extent ' + extent + ' -noise - dwi.mif tmp.mif | mrcalc - -finite - 0 -if lowbnoisemap.mif')
        run.function(os.remove,'tmp.mif')
        run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwirc.mif working.mif')
elif app.args.rician_lowsnr:
    print("...Beginning Rician correction")
    bvalu = np.unique(np.around(bval, decimals=-1))
    lowbval = [i for i in bvalu if i<=np.int(app.args.rician_lowsnr)]
    lowbvalstr = ','.join(str(i) for i in lowbval)
    run.command('dwiextract -shell ' + lowbvalstr + ' dwi.mif dwilowb.mif')
    run.command('dwidenoise -extent ' + extent + ' -noise - dwilowb.mif tmp.mif | mrcalc - -finite - 0 -if lowbnoisemap.mif')
    run.function(os.remove,'tmp.mif')
    run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwirc.mif working.mif')

# b0 normalisation
if app.args.normalise:
    print("...Beginning normalisation")
    if len(DWInlist) == 1:
        run.command('dwinormalise working.mif CSFmask.nii dwinm.mif')
    else:
        miflist = []
        for idx,i in enumerate(DWInlist):
            run.command('mrconvert -coord 3 ' + idxlist[idx] + ' working.mif dwiprenm' + str(idx) + '.mif')
            run.command('dwinormalise dwiprenm' + str(idx) + '.mif CSFmask.nii dwinm' + str(idx) + '.mif')
            miflist.append('dwinm' + str(idx) + '.mif')
            DWImif = ' '.join(miflist)
        run.command('mrcat -axis 3 ' + DWImif + ' dwinm.mif')
    run.function(os.remove,'working.mif')
    run.command('mrconvert dwinm.mif working.mif')

run.command('mrconvert -force -datatype float32le working.mif dwi_designer.nii')
run.function(os.remove,'working.mif')

if app.args.DTIparams or app.args.DKIparams or app.args.WMTIparams:
    if not os.path.exists(path.fromUser(app.args.output, True)):
        os.makedirs(path.fromUser(app.args.output, True))
    shutil.copyfile(path.toTemp('dwi_designer.nii',True),path.fromUser(app.args.output + '/dwi_designer.nii', True))
    shutil.copyfile(path.toTemp('dwi_designer.bvec',True),path.fromUser(app.args.output + '/dwi_designer.bvec', True))
    shutil.copyfile(path.toTemp('dwi_designer.bval',True),path.fromUser(app.args.output + '/dwi_designer.bval', True))
    
    print("...Beginning tensor estimation")
    os.chdir(designer_root)
    eng = matlab.engine.start_matlab()

    outliers=0
    if app.args.outliers:
        outliers=1
    DTI=0
    if app.args.DTIparams:
        DTI=1
    DKI=0
    if app.args.DKIparams:
        DKI=1
    WMTI=0
    if app.args.WMTIparams:
        WMTI=1
    constraints=0
    if app.args.fit_constraints:
        constraints=app.args.fit_constraints
    AKC=0
    if app.args.akc:
        AKC=1
    KCUM=0
    if app.args.kcumulants:
        KCUM=1

    eng.tensorfitting(path.toTemp('',True),path.fromUser(app.args.output, True),outliers,KCUM,DTI,DKI,WMTI,constraints,AKC,DKI_root,nargout=0)
    eng.quit()
    app.gotoTempDir()
else:
    if app.args.datatype:
        run.command('mrconvert -datatype ' + app.args.datatype + ' ' + path.toTemp('dwi_designer.nii',True) + ' ' + path.fromUser(app.args.output + '.nii',True))
    else:
        shutil.copyfile(path.toTemp('dwi_designer.nii',True),path.fromUser(app.args.output + '.nii',True))
    shutil.copyfile(path.toTemp('dwi_designer.bvec',True),path.fromUser(app.args.output + '.bvec',True))
    shutil.copyfile(path.toTemp('dwi_designer.bval',True),path.fromUser(app.args.output + '.bval',True))


app.complete()

