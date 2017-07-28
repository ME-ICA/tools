#!/usr/bin/env python
__version__="beta1"
welcome_block="""
# AFNI HCP Preprocessing, Version %s
#
#
# topup.py (c) 2014 Prantik Kundu
# Use FSL TOPUP to estimate field map
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
fmopts=OptionParser()
fmopts.add_option('-d',"--phase_dir",dest='phase_dir',help="Main phase encoding direction x=1 -x=-1 y=2 -y=-2. ex: --phase_dir=-1 ",default=None)
fmopts.add_option('-A',"--smphase_sample",dest='smphase_sample',help="Sample of fMRI data with main phase encoding",default=None)
fmopts.add_option('-B',"--opphase_sample",dest='opphase_sample',help="Sample of fMRI data with opposite phase encoding",default=None)
fmopts.add_option('-r',"--readout_time",dest='readout_time',help="Readout time",default=None)
fmopts.add_option('-s',"--prefix",dest='setname',help="Name of this fieldmap",default=None)

runopts=OptionGroup(fmopts,"Run options")
runopts.add_option('',"--cpus",dest='cpus',help="Maximum number of CPUs (OpenMP threads) to use. Default 2.",default='2')
runopts.add_option('',"--test_proc",action="store_true",dest='test_proc',help="Align and preprocess 1 dataset then exit, for testing",default=False)
runopts.add_option('',"--script_only",action="store_true",dest='script_only',help="Generate script only, then exit",default=0)
runopts.add_option('',"--OVERWRITE",dest='overwrite',action="store_true",help="If hcpp.xyz directory exists, overwrite. ",default=False)
fmopts.add_option_group(runopts)

(options,args) = fmopts.parse_args()

#Prepare script
startdir=rstrip(popen('pwd').readlines()[0])
hcppdir=os.path.dirname(sys.argv[0])
sl = []	#Script command list
sl.append('#'+" ".join(sys.argv).replace('"',r"\""))
sl.append(welcome_block)
osf='.nii.gz' #Using NIFTI outputs
setname=options.setname

#Parse field map options
if options.readout_time and options.phase_dir and options.opphase_sample:
	make_fm = True
	fm_correct = True
	print "++ Will compute field map from oppositely encoded data"
else:
	print "*+ Field map options specified, but could not interpret! Ask for help."
	sys.exit()

#Detect AFNI direcotry
afnidir = os.path.dirname(os.popen('which 3dSkullStrip').readlines()[0])

#Prepare script and enter hcpp directory
logcomment("Set up script run environment",level=1)
#sl.append('module load afni fsl gsl')
sl.append('set -e')
sl.append('export OMP_NUM_THREADS=%s' % (options.cpus))
sl.append('export MKL_NUM_THREADS=%s' % (options.cpus))
sl.append('export DYLD_FALLBACK_LIBRARY_PATH=%s' % (afnidir))
sl.append('export AFNI_3dDespike_NEW=YES')

if options.overwrite: 
	sl.append('rm -rf hcpp-fm.%s' % (setname))
else: 
	sl.append("if [[ -e hcpp-fm.%s ]]; then echo HCPP-FM directory exists, exiting; exit; fi" % (setname))
sl.append('mkdir -p hcpp-fm.%s' % (setname))
sl.append("cp _hcpp-fm_%s.sh hcpp-fm.%s/" % (setname,setname))
sl.append("cd hcpp-fm.%s" % setname)
thecwd= "%s/hcpp-fm.%s" % (getcwd(),setname)

logcomment("Perform distortion correction using FSL TOPUP",level=1)
if fm_correct:
	#corrected_out = '%s_fmc.nii.gz' % (dsprefix(datasets_in[0]))
	readout_time = float(options.readout_time)
	fmarray=[[0,0,0,readout_time],[0,0,0,readout_time]]
	fm_axis = abs(int(options.phase_dir))-1
	fm_dir = int(options.phase_dir)/abs(int(options.phase_dir))
	fmarray[0][fm_axis]=fm_dir
	fmarray[1][fm_axis]=fm_dir*-1
	#Save acqparams.txt
	sl.append("printf -- '%s\\n' > acqparams.txt" % '\\n'.join(['\\t'.join([str(vv) for vv in ll]) for ll in fmarray ]))
	#Get TOPUP field map
	sl.append('fslroi %s/%s ./fm_smphs.nii.gz 0 1' % (startdir,options.smphase_sample) )
	sl.append("fslroi %s/%s ./fm_opphs.nii.gz 0 1" % (startdir,options.opphase_sample) )
	sl.append('fslmerge -t fm_both fm_smphs.nii.gz fm_opphs.nii.gz' )
	sl.append('topup --imain=fm_both --datain=acqparams.txt --config=b02b0.cnf --out=fm')
	for ii in range(3): sl.append('echo \'%s\' >> fm_movpar.txt' % '  '.join(['0']*6) )
	#sl.append('applytopup --imain=%s/%s --datain=acqparams.txt --topup=fm --out=%s/%s --inindex=1 --method=jac' % (startdir,datasets_in[0],startdir,corrected_out))

#Write the preproc script and execute it
ofh = open('_hcpp-fm_%s.sh' % setname ,'w')
print "++ Writing script file: _hcpp-fm_%s.sh" % (setname)
ofh.write("\n".join(sl)+"\n")
ofh.close()
if not options.script_only: 
	print "++ Executing script file: _hcpp-fm_%s.sh" % (setname)
	system('bash _hcpp-fm_%s.sh' % setname)
