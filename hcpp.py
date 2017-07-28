#!/usr/bin/env python
__version__="beta1"
welcome_block="""
# AFNI HCP Preprocessing, Version %s
#
# Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S., 
# Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional 
# connectivity mapping using multiecho fMRI. PNAS (2013).
#
# hcpp.py (c) 2014 Prantik Kundu
# PROCEDURE 1 : Preprocess HCP datasets
# -Check arguments, input filenames, and filesystem for dependencies
# -Calculation of motion parameters
# -Application of motion correction and coregistration parameters
# -Misc. EPI preprocessing (temporal alignment, smoothing, etc) in appropriate order
""" % (__version__)

import sys
import commands
from re import split as resplit
import re
from os import system,getcwd,mkdir,chdir,popen
import os.path
from string import rstrip,split
from optparse import OptionParser,OptionGroup,SUPPRESS_HELP

#Filename parser for NIFTI and AFNI files
def dsprefix(idn):
	def prefix(datasetname):
		return split(datasetname,'+')[0]
	if len(split(idn,'.'))!=0:
		if split(idn,'.')[-1]=='HEAD' or split(idn,'.')[-1]=='BRIK' or split(idn,'.')[-2:]==['BRIK','gz']:
			return prefix(idn)
		elif split(idn,'.')[-1]=='nii' and not split(idn,'.')[-1]=='nii.gz':
			return '.'.join(split(idn,'.')[:-1])
		elif split(idn,'.')[-2:]==['nii','gz']:
			return '.'.join(split(idn,'.')[:-2])
		else:
			return prefix(idn)
	else:
		return prefix(idn)

def dssuffix(idna):
	suffix = idna.split(dsprefix(idna))[-1]
	#print suffix
	spl_suffix=suffix.split('.')
	#print spl_suffix
	if len(spl_suffix[0])!=0 and spl_suffix[0][0] == '+': return spl_suffix[0]
	else: return suffix

def getdsname(e_ii,prefixonly=False):
	if shorthand_dsin: dsname = '%s%s%s%s' % (prefix,datasets[e_ii],trailing,isf)
	else: dsname =  datasets_in[e_ii]
	if prefixonly: return dsprefix(dsname)
	else: return dsname

def logcomment(comment,level=3): 
	majmark='\n'
	leading='--------'
	if level==3: majmark=''
	if level==1: 
		leading="+* "
		sl.append("""echo "\n++++++++++++++++++++++++" """)
		majmark=''
	sl.append("""%secho %s"%s" """ % (majmark,leading,comment))

#Configure options and help dialog
parser=OptionParser()
parser.add_option('-d',"",dest='dsinputs',help="ex: -d REST1.nii.gz",default='')
parser.add_option('-a',"",dest='anat',help="ex: -a mprage.nii.gz  Anatomical dataset (optional)",default='')
parser.add_option('-b',"",dest='basetime',help="ex: -b 10s OR -b 10v  Time to steady-state equilibration in seconds(s) or volumes(v). Default 0. ",default='0')
parser.add_option('',"--MNI",dest='mni',action='store_true',help="Warp to MNI space using high-resolution template",default=False)
extopts=OptionGroup(parser,"Additional processing options")
extopts.add_option('',"--qwarp",dest='qwarp',action='store_true',help="Nonlinear anatomical normalization to MNI (or --space template) using 3dQWarp, after affine",default=False)
extopts.add_option('',"--fres",dest='fres',help="Specify functional voxel dim. in mm (iso.) for resampling during preprocessing. Default none. ex: --fres=2.5", default=False)
extopts.add_option('',"--space",dest='space',help="Path to specific standard space template for affine anatomical normalization",default=False)
extopts.add_option('',"--no_skullstrip",action="store_true",dest='no_skullstrip',help="Anatomical is already intensity-normalized and skull-stripped (if -a provided)",default=False)
extopts.add_option('',"--no_despike",action="store_true",dest='no_despike',help="Do not de-spike functional data. Default is to de-spike, recommended.",default=False)
extopts.add_option('',"--no_axialize",action="store_true",dest='no_axialize',help="Do not re-write dataset in axial-first order. Default is to axialize, recommended.",default=False)
extopts.add_option('',"--smooth",dest='FWHM',help="Data FWHM smoothing (3dBlurInMask). Default off. ex: --smooth 3mm ",default='0mm')
extopts.add_option('',"--align_base",dest='align_base',help="Explicitly specify base dataset for volume registration",default='')
extopts.add_option('',"--TR",dest='TR',help="The TR. Default read from input dataset header",default='')
extopts.add_option('',"--tpattern",dest='tpattern',help="Slice timing (i.e. alt+z, see 3dTshift -help). Default from header. (N.B. This is important!)",default='')
extopts.add_option('',"--align_args",dest='align_args',help="Additional arguments to anatomical-functional co-registration routine",default='')
parser.add_option_group(extopts)
fmopts=OptionGroup(parser,"Fieldmap correction options")
fmopts.add_option('',"--fm_refloc",dest='fm_refloc',help="If field map already computed, indicate directory",default=None)
fmopts.add_option('',"--fm_refdir",dest='fm_refdir',help="If field map already computed, indicate its phase encoding, see --phase_dir options",default=None)
fmopts.add_option('',"--phase_dir",dest='phase_dir',help="Phase encoding direction x=1 -x=-1 y=2 -y=-2. ex: --phase_dir=-1 ",default=None)
fmopts.add_option('',"--opphase_sample",dest='opphase_sample',help="Sample of fMRI data with opposite phase encoding",default=None)
fmopts.add_option('',"--readout_time",dest='readout_time',help="Readout time",default=None)
parser.add_option_group(fmopts)
dnopts=OptionGroup(parser,"Basic denoising options")
dnopts.add_option('',"--no_moreg",dest='no_moreg',action='store_true',help="No 0th order motion regression",default=False)
dnopts.add_option('',"--no_modtreg",dest='no_modtreg',action='store_true',help="No 1st order motion regression",default=False)
dnopts.add_option('',"--highpass",dest='highpass',help="High pass filter",default='0')
dnopts.add_option('',"--lowpass",dest='lowpass',help="Low pass filter",default='0')
dnopts.add_option('',"--detrend",dest='detrend',help="Order for polynomial detrend",default='-1')
parser.add_option_group(dnopts)
runopts=OptionGroup(parser,"Run optipns")
runopts.add_option('',"--prefix",dest='prefix',help="Prefix for final HCPP output datasets.",default='')
runopts.add_option('',"--cpus",dest='cpus',help="Maximum number of CPUs (OpenMP threads) to use. Default 2.",default='2')
runopts.add_option('',"--label",dest='label',help="Label to tag HCPP analysis folder.",default='')
runopts.add_option('',"--test_proc",action="store_true",dest='test_proc',help="Align and preprocess 1 dataset then exit, for testing",default=False)
runopts.add_option('',"--script_only",action="store_true",dest='script_only',help="Generate script only, then exit",default=0)
runopts.add_option('',"--keep_int",action="store_true",dest='keep_int',help="Keep preprocessing intermediates. Default delete.",default=False)
runopts.add_option('',"--skip_check",action="store_true",dest='skip_check',help="Skip dependency checks during initialization.",default=False)
runopts.add_option('',"--FINALIZE",dest='finalize',action="store_true",help="If hcpp already run up to VR, then continue.", default=False)
runopts.add_option('',"--OVERWRITE",dest='overwrite',action="store_true",help="If hcpp.xyz directory exists, overwrite. ",default=False)
parser.add_option_group(runopts)
(options,args) = parser.parse_args()

#Welcome line
print """\n-- HCP preprocessing based on ME-ICA Pipeline %s --

Please cite: 
Kundu, P., Brenowitz, N.D., Voon, V., Worbe, Y., Vertes, P.E., Inati, S.J., Saad, Z.S., 
Bandettini, P.A. & Bullmore, E.T. Integrated strategy for improving functional 
connectivity mapping using multiecho fMRI. PNAS (2013).

Kundu, P., Inati, S.J., Evans, J.W., Luh, W.M. & Bandettini, P.A. Differentiating 
BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage (2011).
""" % (__version__)

#Parse dataset input names
if options.dsinputs=='' or options.TR==0:
	#dep_check()
	print "*+ Need at least dataset inputs and TE. Try hcpp.py -h"
	sys.exit()
if os.path.abspath(os.path.curdir).__contains__('hcpp.'):
	print "*+ You are inside a HCPP directory! Please leave this directory and rerun."
	sys.exit()

"""Hacks to get this into single-echo only mode for HCP PP"""
options.coreg_mode='aea'
options.mask_mode='func'
shorthand_dsin = False
datasets_in=[options.dsinputs]
datasets = ['1']
prefix = dsprefix(datasets_in[0])
isf = dssuffix(datasets_in[0])
if '.nii' in isf: isf='.nii'
trailing=''
setname=prefix+options.label
tes=['default']

#Prepare script
startdir=rstrip(popen('pwd').readlines()[0])
hcppdir=os.path.dirname(sys.argv[0])
headsl = []
sl = []	#Script command list
headsl.append('#'+" ".join(sys.argv).replace('"',r"\""))
headsl.append(welcome_block)
osf='.nii.gz' #Using NIFTI outputs

#Check if input files exist
notfound=0
for ds_ii in range(len(datasets)): 
	if commands.getstatusoutput('3dinfo %s' % (getdsname(ds_ii)))[0]!=0:
		print "*+ Can't find/load dataset %s !" % (getdsname(ds_ii))
		notfound+=1
if options.anat!='' and commands.getstatusoutput('3dinfo %s' % (options.anat))[0]!=0:
	print "*+ Can't find/load anatomical dataset %s !" % (options.anat)
	notfound+=1
if notfound!=0:
	print "++ EXITING. Check dataset names."
	sys.exit()

#Parse timing arguments
if options.TR!='':tr=float(options.TR)
else: 
	tr=float(os.popen('3dinfo -tr %s' % (getdsname(0))).readlines()[0].strip())
	options.TR=str(tr)
if 'v' in str(options.basetime): 
	basebrik = int(options.basetime.strip('v'))
else:
	timetoclip=0
	timetoclip = float(options.basetime.strip('s'))
	basebrik=int(round(timetoclip/tr))

#Misc. command parsing
if options.mni: options.space='MNI_caez_N27+tlrc'
if options.qwarp and (options.anat=='' or not options.space):
	print "*+ Can't specify Qwarp nonlinear coregistration without anatomical and SPACE template!"
	sys.exit()

if not options.mask_mode in ['func','anat','template']:
	print "*+ Mask mode option '%s' is not recognized!" % options.mask_mode
	sys.exit()

#Parse field map options
if [None]==list(set([options.fm_refloc,options.fm_refdir, options.readout_time,options.phase_dir,options.opphase_sample])):
	make_fm = False
	fm_correct = False
elif options.fm_refloc or options.fm_refdir or options.readout_time or options.phase_dir or options.opphase_sample:
	if options.fm_refloc and options.fm_refdir and options.phase_dir:
		make_fm = False
		fm_correct = True
		print "++ Using reference field map from %s" % options.fm_refloc
	elif options.readout_time and options.phase_dir and options.opphase_sample:
		make_fm = True
		fm_correct = True
		print "++ Will compute field map from oppositely encoded data"
else:
	print "*+ Field map options specified, but could not interpret! Ask for help."
	sys.exit()

#Parse alignment options
if options.coreg_mode == 'aea': options.t2salign=False
elif 'lp' in options.coreg_mode : options.t2salign=True
align_base = basebrik
align_interp='cubic'
oblique_epi_read = 0 
oblique_anat_read = 0
zeropad_opts = " -I %s -S %s -A %s -P %s -L %s -R %s " % (tuple([1]*6))
if options.anat!='':
	oblique_anat_read = int(os.popen('3dinfo -is_oblique %s' % (options.anat)).readlines()[0].strip())
	epicm = [float(coord) for coord in os.popen("3dCM %s" % (getdsname(0))).readlines()[0].strip().split()]
	anatcm = [float(coord) for coord in os.popen("3dCM %s" % (options.anat)).readlines()[0].strip().split()]
	maxvoxsz = float(os.popen("3dinfo -dk %s" % (getdsname(0))).readlines()[0].strip())
	deltas = [abs(epicm[0]-anatcm[0]),abs(epicm[1]-anatcm[1]),abs(epicm[2]-anatcm[2])]
	cmdist = 20+sum([dd**2. for dd in deltas])**.5
	cmdif =  max(abs(epicm[0]-anatcm[0]),abs(epicm[1]-anatcm[1]),abs(epicm[2]-anatcm[2]))
	addslabs = abs(int(cmdif/maxvoxsz))+10
   	zeropad_opts=" -I %s -S %s -A %s -P %s -L %s -R %s " % (tuple([addslabs]*6))
oblique_epi_read = int(os.popen('3dinfo -is_oblique %s' % (getdsname(0))).readlines()[0].strip())
if oblique_epi_read or oblique_anat_read: 
	oblique_mode = True
	sl.append("echo Oblique data detected.")
else: oblique_mode = False
if options.fres:
	if options.qwarp: qwfres="-dxyz %s" % options.fres
	else: alfres = "-mast_dxyz %s" % options.fres
else: 
	if options.qwarp: qwfres="-dxyz ${voxsize}" #See section called "Preparing functional masking for this ME-EPI run"
	else: alfres="-mast_dxyz ${voxsize}"
if options.anat=='' and options.mask_mode!='func':
	print "*+ Can't do anatomical-based functional masking without an anatomical!"
	sys.exit()

#Detect AFNI direcotry
afnidir = os.path.dirname(os.popen('which 3dSkullStrip').readlines()[0])

#Detect finalize mode
finalize_mode = False	
if options.finalize and os.path.exists("hcpp.%s/e1_vr.nii.gz" % setname) and os.path.exists("hcpp.%s/eBvrmask.nii.gz" % setname) and os.path.exists("hcpp.%s/gms.1D" % setname) and os.path.exists("hcpp.%s/finalize_ready" % setname)  :
	finalize_mode = True
	options.overwrite = False
	print "Doing final steps only!"

#Prepare script and enter hcpp directory
logcomment("Set up script run environment",level=1)
#headsl.append('module unload gcc; module load afni fsl/5.0.7 gsl')
headsl.append('set -e')
headsl.append('export OMP_NUM_THREADS=%s' % (options.cpus))
headsl.append('export MKL_NUM_THREADS=%s' % (options.cpus))
headsl.append('export DYLD_FALLBACK_LIBRARY_PATH=%s' % (afnidir))
headsl.append('export AFNI_3dDespike_NEW=YES')
headsl.append('module load afni gsl')

if options.overwrite: 
	headsl.append('rm -rf hcpp.%s' % (setname))
elif not finalize_mode:
	headsl.append("if [[ -e hcpp.%s ]]; then echo HCPP directory exists, exiting; exit; fi" % (setname))
headsl.append('mkdir -p hcpp.%s' % (setname))
if finalize_mode: scriptname = "_hcpp_%s_finalize.sh" % (setname)
else: scriptname = "_hcpp_%s.sh" % (setname)
headsl.append("cp %s hcpp.%s/" % (scriptname,setname))
headsl.append("cd hcpp.%s" % setname)
thecwd= "%s/hcpp.%s" % (getcwd(),setname)

ica_datasets = sorted(datasets)

#Parse anatomical processing options, process anatomical
if options.anat != '':
	logcomment("Deoblique, unifize, skullstrip, and/or autobox anatomical, in starting directory (may take a little while)", level=1)
	nsmprage = options.anat
	anatprefix=dsprefix(nsmprage)
	pathanatprefix="%s/%s" % (startdir,anatprefix)
	if oblique_mode:
		sl.append("if [ ! -e %s_do.nii.gz ]; then 3dWarp -overwrite -prefix %s_do.nii.gz -deoblique %s/%s; fi" % (pathanatprefix,pathanatprefix,startdir,nsmprage))
		nsmprage="%s_do.nii.gz" % (anatprefix)
	if not options.no_skullstrip: 
		sl.append("if [ ! -e %s_ns.nii.gz ]; then 3dUnifize -overwrite -prefix %s_u.nii.gz %s/%s; medline=`3dBrickStat -percentile 99 1 99 %s_do.nii.gz`; medlinea=($medline); p50=${medlinea[1]}; 3dcalc -a %s_do.nii.gz -overwrite -b %s_u.nii.gz -expr \"b*step(a-$p50/20)\" -prefix %s_uc.nii.gz; 3dSkullStrip -overwrite -prefix %s_ns.nii.gz -input %s_uc.nii.gz; 3dAutobox -overwrite -prefix %s_ns.nii.gz %s_ns.nii.gz; fi" % (pathanatprefix,pathanatprefix,startdir,nsmprage,pathanatprefix,pathanatprefix,pathanatprefix,pathanatprefix,pathanatprefix,pathanatprefix,pathanatprefix, pathanatprefix))
		nsmprage="%s_ns.nii.gz" % (anatprefix)

sl.append('module unload gcc gsl python fsl afni; module load fsl/5.0.7')

logcomment("Perform distortion correction using FSL TOPUP",level=1)
if fm_correct:
	corrected_out = '%s_fmc.nii.gz' % (dsprefix(datasets_in[0]))
	if make_fm:
		readout_time = float(options.readout_time)
		fmarray=[[0,0,0,readout_time],[0,0,0,readout_time]]
		fm_axis = abs(int(options.phase_dir))-1
		fm_dir = int(options.phase_dir)/abs(int(options.phase_dir))
		fmarray[0][fm_axis]=fm_dir
		fmarray[1][fm_axis]=fm_dir*-1
		#Save acqparams.txt
		sl.append("printf -- '%s\\n' > acqparams.txt" % '\\n'.join(['\\t'.join([str(vv) for vv in ll]) for ll in fmarray ]))
		#Get TOPUP field map
		sl.append('fslroi %s/%s ./fm_smphs.nii.gz 0 1' % (startdir,datasets_in[0]) )
		sl.append("fslroi %s/%s ./fm_opphs.nii.gz 0 1" % (startdir,options.opphase_sample) )
		sl.append('fslmerge -t fm_both fm_smphs.nii.gz fm_opphs.nii.gz' )
		sl.append('topup --imain=fm_both --datain=acqparams.txt --config=b02b0.cnf --out=fm' )
		sl.append('applytopup --imain=%s/%s --datain=acqparams.txt --topup=fm --out=%s/%s --inindex=1 --method=jac' % (startdir,datasets_in[0],startdir,corrected_out))
	else: #Assume ref fm indicated
		if options.fm_refloc[0]=='/': fm_refloc = options.fm_refloc
		else: fm_refloc = '%s/%s' % (startdir,options.fm_refloc)
		if int(options.fm_refdir)==int(options.phase_dir): 
			fm_coefloc = fm_refloc
			fm_acqloc = fm_refloc
			fm_refindex=1
		elif int(options.fm_refdir)==-1*int(options.phase_dir): 
			fm_coefloc = fm_refloc
			fm_acqloc = fm_refloc
			fm_refindex=2
		elif abs(int(options.fm_refdir))!=abs(int(options.phase_dir)):
			#existing_fm = open('%s/acqparams.txt' % options.fm_refloc).readlines()
			sl.append("cp %s/acqparams.txt . " %  fm_refloc )
			readout_time = float(options.readout_time)
			fmline=[0,0,0,readout_time]
			fm_axis = abs(int(options.phase_dir))-1
			fm_dir = int(options.phase_dir)/abs(int(options.phase_dir))
			fmline[fm_axis]=fm_dir
			sl.append("printf -- '%s\\n' >> acqparams.txt" % '\\t'.join([str(vv) for vv in fmline]) )
			fm_coefloc = fm_refloc
			fm_acqloc = '.'
			fm_refindex=3
		sl.append('applytopup --imain=%s/%s --datain=%s/acqparams.txt --topup=%s/fm --out=%s/%s --inindex=%i --method=jac' % (startdir,datasets_in[0],fm_acqloc,fm_coefloc,startdir,corrected_out,fm_refindex))
	datasets_in[0] = corrected_out

sl.append('module load afni gsl')

# Copy in functional datasets as NIFTI (if not in NIFTI already), calculate rigid body alignment
vrbase=getdsname(0,True)
logcomment("Copy in functional datasets, reset NIFTI tags as needed", level=1)
for e_ii in range(len(datasets)):
	ds = datasets[e_ii]
	sl.append("3dcalc -a %s/%s -expr 'a' -prefix ./%s.nii" % (startdir,getdsname(e_ii),getdsname(e_ii,True) )   )
	if '.nii' in isf: 
		sl.append("nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles ./%s.nii -overwrite" % (  getdsname(e_ii,True)  ))
isf = '.nii'

logcomment("Calculate and save motion and obliquity parameters, despiking first if not disabled, and separately save and mask the base volume",level=1)
#Determine input to volume registration
vrAinput = "./%s%s" % (vrbase,isf)
#Compute obliquity matrix
if oblique_mode: 
	if options.anat!='': 
		if not options.no_axialize:
			sl.append("3daxialize -prefix ob_sam.nii.gz %s[0]" % vrAinput)
			ob_sam = 'ob_sam.nii.gz'
		else:
			ob_sam = "%s[0]" % vrAinput
		sl.append("3dWarp -verb -card2oblique %s -overwrite  -newgrid 1.000000 -prefix ./%s_ob.nii.gz %s/%s | \grep  -A 4 '# mat44 Obliquity Transformation ::'  > %s_obla2e_mat.1D" % (ob_sam,anatprefix,startdir,nsmprage,prefix))
	else: sl.append("3dWarp -overwrite -prefix %s -deoblique %s" % (vrAinput,vrAinput))
#Despike and axialize
if not options.no_despike:
	sl.append("3dDespike -overwrite -prefix ./%s_vrA%s %s "  % (vrbase,osf,vrAinput))
	vrAinput = "./%s_vrA%s" % (vrbase,osf)
if not options.no_axialize: 
	sl.append("3daxialize -overwrite -prefix ./%s_vrA%s %s" % (vrbase,osf,vrAinput))
	vrAinput = "./%s_vrA%s" % (vrbase,osf)
#Set eBbase
external_eBbase=False
if options.align_base!='':
	if options.align_base.isdigit():
		basevol = '%s[%s]' % (vrAinput,options.align_base)
	else:
		basevol = options.align_base
		external_eBbase=True
else: 
	basevol = '%s[%s]' % (vrAinput,basebrik)
sl.append("3dcalc -a %s  -expr 'a' -prefix eBbase.nii.gz "  % (basevol) )
if external_eBbase:
	if oblique_mode: sl.append("3dWarp -overwrite -deoblique eBbase.nii.gz eBbase.nii.gz")
	if not options.no_axialize: sl.append("3daxialize -overwrite -prefix eBbase.nii.gz eBbase.nii.gz")
#Compute motion parameters
sl.append("3dvolreg -overwrite -tshift -quintic  -prefix ./%s_vrA%s -base eBbase.nii.gz -dfile ./%s_vrA.1D -1Dmatrix_save ./%s_vrmat.aff12.1D %s" % \
		  (vrbase,osf,vrbase,prefix,vrAinput))
vrAinput = "./%s_vrA%s" % (vrbase,osf)
sl.append("1dcat './%s_vrA.1D[1..6]{%s..$}' > motion.1D " % (vrbase,basebrik))
e2dsin = prefix+datasets[0]+trailing

logcomment("Preliminary preprocessing of functional datasets: despike, tshift, deoblique, and/or axialize",level=1)
#Do preliminary preproc for this run
if shorthand_dsin: datasets.sort()
for echo_ii in range(len(datasets)):
	#Determine dataset name
	echo = datasets[echo_ii]
	indata = getdsname(echo_ii)
	dsin = 'e'+echo
	if echo_ii==0: e1_dsin = dsin
	logcomment("Preliminary preprocessing dataset %s of TE=%sms to produce %s_ts+orig" % (indata,str(tes[echo_ii]),dsin) )
	#Pre-treat datasets: De-spike, RETROICOR in the future?
	intsname = '%s%s' % (dsprefix(indata),isf)
	if not options.no_despike:
		intsname = "./%s_pt.nii.gz" % dsprefix(indata)
		sl.append("3dDespike -overwrite -prefix %s %s%s" % (intsname,dsprefix(indata),isf))
	#Time shift datasets
	if options.tpattern!='':
		tpat_opt = ' -tpattern %s ' % options.tpattern
	else:
		tpat_opt = ''
	sl.append("3dTshift -heptic %s -prefix ./%s_ts+orig %s" % (tpat_opt,dsin,intsname) )
	if oblique_mode and options.anat=="":
		sl.append("3dWarp -overwrite -deoblique -prefix ./%s_ts+orig ./%s_ts+orig" % (dsin,dsin))
	#Axialize functional dataset
	if not options.no_axialize:
		sl.append("3daxialize  -overwrite -prefix ./%s_ts+orig ./%s_ts+orig" % (dsin,dsin))
	if oblique_mode: sl.append("3drefit -deoblique -TR %s %s_ts+orig" % (options.TR,dsin))
	else: sl.append("3drefit -TR %s %s_ts+orig" % (options.TR,dsin))
#Compute grand mean scaling factor
sl.append("3dBrickStat -mask eBbase.nii.gz -percentile 50 1 50 %s_ts+orig[%i] > gms.1D" % (e1_dsin,basebrik))
sl.append("gms=`cat gms.1D`; gmsa=($gms); p50=${gmsa[1]}")

#Make base volume
sl.append("3dSkullStrip -prefix eBbase_ss.nii.gz -input eBbase.nii.gz")
sl.append("3dcopy eBbase.nii.gz eBbase")

# Calculate affine anatomical warp if anatomical provided, then combine motion correction and coregistration parameters 
if options.anat!='':
	#Copy in anatomical and make sure its in +orig space
	logcomment("Copy anatomical into HCPP directory and process warps",level=1)
	sl.append("cp %s/%s* ." % (startdir,nsmprage))
	abmprage = nsmprage
	refanat = nsmprage
	if options.space:
		sl.append("afnibinloc=`which 3dSkullStrip`")
		if '/' in options.space: 
			sl.append("ll=\"%s\"; templateloc=${ll%%/*}/" % options.space)
			options.space=options.space.split('/')[-1]
  		else:
			sl.append("templateloc=${afnibinloc%/*}")
		atnsmprage = "%s_at.nii.gz" % (dsprefix(nsmprage))
		if not dssuffix(nsmprage).__contains__('nii'): sl.append("3dcalc -float -a %s -expr 'a' -prefix %s.nii.gz" % (nsmprage,dsprefix(nsmprage)))
		logcomment("If can't find affine-warped anatomical, copy native anatomical here, compute warps (takes a while) and save in start dir. ; otherwise link in existing files")
		sl.append("if [ ! -e %s/%s ]; then \@auto_tlrc -no_ss -init_xform AUTO_CENTER -base ${templateloc}/%s -input %s.nii.gz -suffix _at" % (startdir,atnsmprage,options.space,dsprefix(nsmprage)))
		sl.append("cp %s.nii %s" % (dsprefix(atnsmprage),startdir))
		sl.append("gzip -f %s/%s.nii" % (startdir,dsprefix(atnsmprage)))
		sl.append("else ln -s %s/%s ." % (startdir,atnsmprage))
		refanat = '%s/%s' % (startdir,atnsmprage)
		sl.append("fi")
		sl.append("3dcopy %s/%s.nii.gz %s" % (startdir,dsprefix(atnsmprage),dsprefix(atnsmprage)))
		sl.append("3drefit -view orig %s+tlrc " % dsprefix(atnsmprage) )
		sl.append("3dAutobox -prefix ./abtemplate.nii.gz ${templateloc}/%s" % options.space)
		abmprage = 'abtemplate.nii.gz'
		if options.qwarp:
			logcomment("If can't find non-linearly warped anatomical, compute, save back; otherwise link")
			nlatnsmprage="%s_atnl.nii.gz" % (dsprefix(nsmprage))
			sl.append("if [ ! -e %s/%s ]; then " % (startdir,nlatnsmprage))
			logcomment("Compute non-linear warp to standard space using 3dQwarp (get lunch, takes a while) ")
			sl.append("3dUnifize -overwrite -GM -prefix ./%su.nii.gz %s/%s" % (dsprefix(atnsmprage),startdir,atnsmprage))  
			sl.append("3dQwarp -iwarp -overwrite -resample -useweight -blur 2 2 -duplo -workhard -base ${templateloc}/%s -prefix %s/%snl.nii.gz -source ./%su.nii.gz" % (options.space,startdir,dsprefix(atnsmprage),dsprefix(atnsmprage)))
			sl.append("fi")
			sl.append("ln -s %s/%s ." % (startdir,nlatnsmprage))
			refanat = '%s/%snl.nii.gz' % (startdir,dsprefix(atnsmprage))
	
	#Set anatomical reference for anatomical-functional co-registration
	if oblique_mode: alnsmprage = "./%s_ob.nii.gz" % (anatprefix)
	else: alnsmprage = "%s/%s" % (startdir,nsmprage)
	if options.coreg_mode=='lp-t2s': 
		logcomment("Using alignp_mepi_anat.py to drive T2*-map weighted anatomical-functional coregistration")
		ama_alnsmprage = alnsmprage
		if not options.no_axialize:
			ama_alnsmprage = os.path.basename(alnsmprage)
			sl.append("3daxialize -overwrite -prefix ./%s %s" % (ama_alnsmprage,alnsmprage))
		t2salignpath = 'hcpp.libs/alignp_mepi_anat.py'
		sl.append("%s %s -t t2svm_ss.nii.gz -a %s -p mepi %s" % \
			(sys.executable, '/'.join([hcppdir,t2salignpath]),ama_alnsmprage,options.align_args))
		sl.append("cp alignp.mepi/mepi_al_mat.aff12.1D ./%s_al_mat.aff12.1D" % anatprefix)
	elif options.coreg_mode=='aea':
		logcomment("Using AFNI align_epi_anat.py to drive anatomical-functional coregistration ")
		sl.append("3dcopy %s ./ANAT_ns+orig " % alnsmprage)
		#sl.append("align_epi_anat.py -anat2epi -volreg off -tshift off -deoblique off -epi_strip None -anat_has_skull no -save_script aea_anat_to_ocv.tcsh -anat ANAT_ns+orig -epi eBbase+orig -epi_base 0 %s" % (options.align_args) )
		sl.append("3dAllineate -overwrite  -source ANAT_ns+orig. -base eBbase_ss.nii.gz -warp aff -source_automask+4  -prefix ANAT_ns_al -1Dmatrix_save ANAT_ns_al_mat.aff12.1D  -master SOURCE %s" % (options.align_args))
		sl.append("cp ANAT_ns_al_mat.aff12.1D %s_al_mat.aff12.1D" % (anatprefix))
	if options.space: 
		tlrc_opt = "%s/%s::WARP_DATA -I" % (startdir,atnsmprage)
		sl.append("cat_matvec -ONELINE %s > %s/%s_ns2at.aff12.1D" % (tlrc_opt,startdir,anatprefix))
	else: tlrc_opt = ""
	if oblique_mode: oblique_opt = "%s_obla2e_mat.1D" % prefix
	else: oblique_opt = ""
	sl.append("cat_matvec -ONELINE  %s %s %s_al_mat.aff12.1D -I > %s_wmat.aff12.1D" % (tlrc_opt,oblique_opt,anatprefix,prefix))
	sl.append("cat_matvec -ONELINE  %s %s %s_al_mat.aff12.1D -I  %s_vrmat.aff12.1D  > %s_vrwmat.aff12.1D" % (tlrc_opt,oblique_opt,anatprefix,prefix,prefix))

else: sl.append("cp %s_vrmat.aff12.1D %s_vrwmat.aff12.1D" % (prefix,prefix))

#Preprocess datasets
if shorthand_dsin: datasets.sort()
logcomment("Extended preprocessing of functional datasets",level=1)
for echo_ii in range(len(datasets)):

	#Determine dataset name
	echo = datasets[echo_ii]
	indata = getdsname(echo_ii)
	dsin = 'e'+echo

	if echo_ii == 0: 
		logcomment("Preparing functional masking for this HCP EPI run",2 )
		if options.anat: almaster="-master %s" % abmprage
		else: almaster=""
		sl.append("3dZeropad %s -prefix eBvrmask.nii.gz eBbase_ss.nii.gz[0]" % (zeropad_opts))
		sl.append("voxsize=`ccalc $(3dinfo -voxvol eBvrmask.nii.gz)**.33`") #Set voxel size
		#Create base mask
		if options.anat and options.space and options.qwarp: 
			sl.append("3dNwarpApply -overwrite -nwarp '%s/%s_WARP.nii.gz' -affter '%s_wmat.aff12.1D' %s %s -source eBvrmask.nii.gz -interp %s -prefix ./eBvrmask.nii.gz " % \
			(startdir,dsprefix(nlatnsmprage),prefix,almaster,qwfres,'NN'))
		elif options.anat:
			sl.append("3dAllineate -overwrite -final %s -%s -float -1Dmatrix_apply %s_wmat.aff12.1D -base eBvrmask.nii.gz -input eBvrmask.nii.gz -prefix ./eBvrmask.nii.gz %s %s" % \
			('NN','NN',prefix,almaster,alfres))

		if options.anat=='':
			logcomment("Trim empty space off of mask dataset and/or resample")
			sl.append("3dAutobox -overwrite -prefix eBvrmask%s eBvrmask%s" % (osf,osf) )
			if options.fres: 
				resstring = "-dxyz %s %s %s" % (options.fres,options.fres,options.fres)
				sl.append("3dresample -overwrite -master eBvrmask.nii.gz %s -input eBvrmask.nii.gz -prefix eBvrmask.nii.gz" % (resstring))
		
		sl.append("3dcalc -float -a eBvrmask.nii.gz -expr 'notzero(a)' -overwrite -prefix eBvrmask.nii.gz")

	#logcomment("Extended preprocessing dataset %s of TE=%sms to produce %s_in.nii.gz" % (indata,str(tes[echo_ii]),dsin),level=2 )
	logcomment("Apply combined normalization/co-registration/motion correction parameter set to %s_ts+orig" % dsin)
	if options.qwarp: sl.append("3dNwarpApply -nwarp '%s/%s_WARP.nii.gz' -affter %s_vrwmat.aff12.1D -master eBvrmask.nii.gz -source %s_ts+orig -interp %s -prefix ./%s_vr%s " % \
			(startdir,dsprefix(nlatnsmprage),prefix,dsin,align_interp,dsin,osf))
	else: sl.append("3dAllineate -final %s -%s -float -1Dmatrix_apply %s_vrwmat.aff12.1D -base eBvrmask%s -input  %s_ts+orig -prefix ./%s_vr%s" % \
		(align_interp,align_interp,prefix,osf,dsin,dsin,osf))
	if echo_ii == 0:
		sl.append("3dTstat -min -prefix ./%s_vr_min%s ./%s_vr%s" % (dsin,osf,dsin,osf) )
		sl.append("3dcalc -a eBvrmask.nii.gz -b %s_vr_min%s -expr 'step(a)*step(b)' -overwrite -prefix eBvrmask.nii.gz " % (dsin,osf))
	
#For finalize mode, need gms.1D, eBvrmask, vr
#Easy way to stage/continue - clear command buffer

	if finalize_mode: 
		sl = []
		sl.append("gms=`cat gms.1D`; gmsa=($gms); p50=${gmsa[1]}")

	if options.FWHM=='0mm': 
		sl.append("3dcalc -float -overwrite -a eBvrmask.nii.gz -b ./%s_vr%s[%i..$] -expr 'step(a)*b' -prefix ./%s_sm%s " % (dsin,osf,basebrik,dsin,osf))
	else: 
		sl.append("3dBlurInMask -overwrite -fwhm %s -mask eBvrmask%s -prefix ./%s_sm%s ./%s_vr%s[%i..$]" % (options.FWHM,osf,dsin,osf,dsin,osf,basebrik))
	sl.append("3dcalc -overwrite -float -overwrite -a ./%s_sm%s -expr \"a*10000/${p50}\" -prefix ./%s_sm%s" % (dsin,osf,dsin,osf))
	sl.append("3dTstat -overwrite -prefix ./%s_mean%s ./%s_sm%s" % (dsin,osf,dsin,osf))

	#Denoising section
	sl.append('1d_tool.py -overwrite -derivative -infile motion.1D -write motion_deriv.1D')
	motocat = []
	if not options.no_moreg: motocat.append('motion.1D')
	if not options.no_modtreg: motocat.append('motion_deriv.1D')
	if len(motocat)>0: 
		sl.append('1dcat %s > motion_regr.1D' % ' '.join(motocat))
	if float(options.detrend)>=0: sl.append("3dDetrend -polort %s -overwrite -prefix ./%s_sm%s ./%s_sm%s " % (options.detrend,dsin,osf,dsin,osf) )
	if float(options.highpass)>=0 and float(options.lowpass)>0: 
		sl.append("3dBandpass -ort motion_regr.1D -prefix ./%s_in.nii %f %f ./%s_sm%s " % (dsin,float(options.highpass),float(options.lowpass),dsin,osf) )
	else: 
		sl.append("3dcalc -a %s_sm%s -expr 'a' -prefix %s_in.nii" % (dsin,osf,dsin))
	
	sl.append("3dTstat -stdev -prefix ./%s_std%s ./%s_in.nii" % (dsin,osf,dsin))
	if options.test_proc: sl.append("exit")
	#if not (options.test_proc or options.keep_int): sl.append("rm -f %s_ts+orig* %s_vr%s %s_sm%s" % (dsin,dsin,osf,dsin,osf))

	#Copy back out
	logcomment("Copying results to start directory",level=1)
	tedflag=''
	if options.prefix!='': prefixedset = '_'.join([setname,options.prefix])
	else: prefixedset = setname
	if options.anat!='' and options.space!=False:
		sl.append("nifti_tool -mod_hdr -mod_field sform_code 2 -mod_field qform_code 2 -infiles %s_in.nii -overwrite" % dsin)
	else:
		sl.append("nifti_tool -mod_hdr -mod_field sform_code 1 -mod_field qform_code 1 -infiles %s_in.nii -overwrite" % dsin)
	sl.append("gzip -f %s_in.nii ; mv %s_in.nii.gz %s/%s_pp.nii.gz" % (dsin,dsin,startdir,prefixedset))

sl.append('cp motion.1D %s/%s_motion.1D' % (startdir,prefixedset))
sl.append('cp motion_deriv.1D %s/%s_motiondt.1D' % (startdir,prefixedset))

#Flush temp files
if not options.keep_int:
	sl.append('mkdir -p _keep')
	sl.append('mv *.1D eBvrmask.nii.gz %s_vr%s  _keep/' % (dsin,osf)  )
	sl.append('rm -f *; mv _keep/* .')
	sl.append('touch finalize_ready')

#Finalize command buffer
sl = headsl + sl

#Write the preproc script and execute it
ofh = open(scriptname ,'w')
print "++ Writing script file: %s" % (scriptname)
ofh.write("\n".join(sl)+"\n")
ofh.close()
if not options.script_only: 
	print "++ Executing script file: %s" % (scriptname)
	system('bash %s' % scriptname)
