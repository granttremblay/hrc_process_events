#! /bin/sh

# 30 January 2002

# This is the official template for pipetool regression test scripts.
# In addition to supporting the "SHORTTEST" option, this script also
# allows the user to run individual subtests from the command line.
# The script will accept a series of test identifiers, generally of
# the form "test1" "test2" ... which are to be run.

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!3, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!4 below.


# !!1
# hrc_process_events.t
# test script for hrc_process_events


# !!2
# syntax:
# hrc_process_events.t [<testid> ... ]
 



######################################################################
# subroutine
# error_exit <message>
# Fatal error exit

error_exit()
{
  echo ""   | tee -a $LOGFILE
  echo "$1" | tee -a $LOGFILE
  echo ""   | tee -a $LOGFILE
  echo "${toolname} : FAIL" | tee -a $LOGFILE
  exit 1
}

######################################################################
# subroutine
# keyfilter infile outfile
# filters out CHECKSUM, Dataset, CREATOR, HISTORY, DATASUM, 
#             ASCDSVER, HISTNUM, and DATE
# To filter additional keywords, add s/KEYWORD/Dataset/g; for each.

keyfilter()
{
  cat $1 | sed -e 's/CHECKSUM/Dataset/g;s/COMMENT/Dataset/g;
  s/DATE/Dataset/g;s/CREATOR/Dataset/g;
  s/HISTORY/Dataset/g; s/DATASUM/Dataset/g;s/ASCDSVER/Dataset/g;
  s/BPIXFILE /Dataset /g;s/DEGAP /Dataset /g;
  s/HYPFILE /Dataset /g;s/TAPRING /Dataset /g;
  s/FLATFILE /Dataset /g;s/SATFILE /Dataset /g;
  s/AMPSFFIL /Dataset /g;
  s/COMMENT /Dataset /g;
  s/HISTNUM/Dataset/g' | \
  grep -v Dataset > $2
  zerotest $2
}

######################################################################
# subroutine
# find_tool <toolname>
# checks that tool exists and is runnable

find_tool()
{
  s1=`type $1`
  s2=`echo $s1 | awk -F" " '{ print $3}'`
  if test -x $s2 ; then
    :
  else
    error_exit "tool $1 not found"
  fi
}

######################################################################
# subroutine
# zerotest <file> 
# Makes sure that file is not 0 length.
# Use this to protect yourself against empty files  (which will 
# 'diff' without error).  This can happen when the input file to
# cat $infile | do_something >> $outfile
# is missing.  This is used by keyfilter(), above.

zerotest()
{
 if test -s $1 ;
 then
   :
 else
   echo "ERROR: file $1 is of zero length" >> $LOGFILE
   #  Indicate failure, but do not exit.
   mismatch=0
 fi
}

######################################################################
# subroutine
# pset_hrc_I()
# pre-set parameters for hrc_process_events

pset_hpe()
{
pset hrc_process_events clobber=yes
pset hrc_process_events tstart="TSTART"
pset hrc_process_events tstop="TSTOP"
pset hrc_process_events verbose=0
pset hrc_process_events start=coarse
pset hrc_process_events stop=sky
pset hrc_process_events rand_seed=1
pset hrc_process_events rand_pix_size=0.5

pset hrc_process_events obsfile=NONE

pset hrc_process_events eventdef="{d:time,s:crsv,s:crsu,s:amp_sf,s:av1,s:av2,s:av3,s:au1,s:au2,s:au3,l:raw,s:chip,l:tdet,f:det,f:sky,s:pha,s:pi,s:sumamps,s:chip_id,x:status}"
pset hrc_process_events  stdlev1="{d:time,s:crsv,s:crsu,s:amp_sf,s:av1,s:av2,s:av3,s:au1,s:au2,s:au3,l:raw,s:chip,l:tdet,f:det,f:sky,s:pha,s:pi,s:sumamps,s:chip_id,x:status}"
pset hrc_process_events simlev1="{l:tick,i:scifr,i:mjf,s:mnf,s:evtctr,s:crsu,s:crsv,s:au1,s:au2,s:au3,s:av1,s:av2,s:av3,s:tdetx,s:tdety,s:pha,s:vstat,s:estat}"
pset hrc_process_events fltlev1="{d:time,s:crsv,s:crsu,s:amp_sf,s:av1,s:av2,s:av3,s:au1,s:au2,s:au3,s:chipx,s:chipy,l:tdetx,l:tdety,s:detx,s:dety,s:x,s:y,s:pha,s:sumamps,s:chip_id,l:status}"
pset hrc_process_events badeventdef="{d:time,s:crsu,s:crsv,s:au1,s:au2,s:au3,s:av1,s:av2,s:av3,s:pha}"
pset hrc_process_events badlev1="{d:time,s:crsu,s:crsv,s:au1,s:au2,s:au3,s:av1,s:av2,s:av3,s:pha}"
pset hrc_process_events hsilev1="{d:time,s:crsu,s:crsv,s:au1,s:au2,s:au3,s:av1,s:av2,s:av3,s:chipx,s:chipy,s:tdetx,s:tdety,s:x,s:y,l:fpz,s:pha,s:vstat,s:estat}"

pset hrc_process_events grid_ratio=0.5
pset hrc_process_events pha_ratio=0.5
pset hrc_process_events wire_charge=0

} # pset_hpe


pset_hrc_I()
{
#------------------------------------------------
# setting up the hrc_process_events parameters :
#------------------------------------------------
echo "setting up parameter file for hrc_I"
echo ""
pset hrc_process_events badfile="$OUTDIR/lev1_hrci_bad_out.fits"
pset hrc_process_events logfile="$OUTDIR/lev1_hrci_out.log"
pset hrc_process_events time_offset=0
pset hrc_process_events instrume=HRC-I
pset hrc_process_events degapfile="COEFF"
pset hrc_process_events do_ratio=yes
pset hrc_process_events do_amp_sf_cor=yes

pset hrc_process_events badpixfile="$INDIR/badpix_in.fits"
pset hrc_process_events acaofffile="$INDIR/hrcf461_000N001_aoff1.fits"
pset hrc_process_events alignmentfile="$INDIR/hrcf461_000N001_soff1.fits"
pset hrc_process_events gainfile=NONE
pset hrc_process_events ADCfile=NONE
pset hrc_process_events ampsfcorfile="$INDIR/hrciAMPSFAug21.fits"
pset hrc_process_events tapfile="$INDIR/hrciD1999-07-23tapringtestN0001.fits"
pset hrc_process_events hypfile="$INDIR/hrciD1999-07-23fptestN0001.fits"
pset hrc_process_events ampsatfile="$INDIR/hrciD1999-07-23sattestN0001.fits"
pset hrc_process_events evtflatfile="$INDIR/hrciD1999-07-22eftestN0001.fits"

pset hrc_process_events cfu1=1.068
pset hrc_process_events cfu2=0
pset hrc_process_events cfv1=1.045
pset hrc_process_events cfv2=0
pset hrc_process_events amp_gain=75
}  # subroutine : pset_hrc_I 


######################################################################
# subroutine
# pset_hrc_S()
# pre-set parameters for hrc_process_events

pset_hrc_S()
{
echo "setting up parameter file for hrc_S"
echo ""
pset hrc_process_events badfile="${OUTDIR}/lev1_hrcs_bad_out.fits"
pset hrc_process_events logfile="${OUTDIR}/lev1_hrcs_out.log"
pset hrc_process_events time_offset=0
pset hrc_process_events instrume=HRC-S
pset hrc_process_events degapfile="${INDIR}/hrcsD1999-11-08gapN0002.fits"
pset hrc_process_events do_ratio=yes
pset hrc_process_events do_amp_sf_cor=yes

pset hrc_process_events badpixfile="${INDIR}/hrcf01246_000N001_bpix1.fits"
pset hrc_process_events acaofffile="${INDIR}/hrcf01246_000N001_aoff1.fits"
pset hrc_process_events alignmentfile="${INDIR}/hrcf01246_000N001_soff1.fits"
pset hrc_process_events gainfile=NONE
pset hrc_process_events ADCfile=NONE
pset hrc_process_events ampsfcorfile="${INDIR}/hrcsAMPSFAug21.fits"
pset hrc_process_events tapfile="NONE"
pset hrc_process_events hypfile="${INDIR}/hrcsD1999-07-23fptestN0001.fits"
pset hrc_process_events ampsatfile="${INDIR}/hrcsD1999-07-23sattestN0001.fits"
pset hrc_process_events evtflatfile="${INDIR}/hrcsD1999-07-22eftestN0001.fits"


pset hrc_process_events cfu1=1.18
pset hrc_process_events cfu2=-0.16
pset hrc_process_events cfv1=1.11
pset hrc_process_events cfv2=-0.1
pset hrc_process_events amp_gain=75
} #subroutine  pset_hrc_S

######################################################################
# subroutine
# pset_hrcS_no_rangelev()
# pre-set parameters for hrc_process_events

pset_hrcS_no_rangelev()
{
echo "setting up parameter file for hrcS_no_rangelev"
echo ""
pset hrc_process_events badfile="${OUTDIR}/lev1_hrcs_bad_out.fits"
pset hrc_process_events logfile="${OUTDIR}/lev1_hrcs_out.log"
pset hrc_process_events time_offset=0
pset hrc_process_events instrume=HRC-S
pset hrc_process_events degapfile="${INDIR}/hrcsD1999-11-08gapN0002.fits"
pset hrc_process_events do_ratio=yes
pset hrc_process_events do_amp_sf_cor=yes

pset hrc_process_events badpixfile="${INDIR}/hrcf01246_000N001_bpix1.fits"
pset hrc_process_events acaofffile="${INDIR}/hrcf01246_000N001_aoff1.fits"
pset hrc_process_events alignmentfile="${INDIR}/hrcf01246_000N001_soff1.fits"
pset hrc_process_events gainfile=NONE
pset hrc_process_events ADCfile=NONE
pset hrc_process_events ampsfcorfile="${INDIR}/hrcsAMPSFAug21.fits"
pset hrc_process_events tapfile="NONE"
pset hrc_process_events hypfile="${INDIR}/hrcsD1999-07-23fptestN0001.fits"
pset hrc_process_events ampsatfile="${INDIR}/hrcsD1999-07-23sattestN0001.fits"
pset hrc_process_events evtflatfile="${INDIR}/hrcsD1999-07-22eftestN0001.fits"

pset hrc_process_events cfu1=1.18
pset hrc_process_events cfu2=-0.16
pset hrc_process_events cfv1=1.11
pset hrc_process_events cfv2=-0.1
pset hrc_process_events amp_gain=75
}  # subroutine :  pset_hrcS_no_rangelev 

######################################################################
# subroutine
# pset_hrcS_low_rangelev()
# pre-set parameters for hrc_process_events

pset_hrcS_low_rangelev()
{
echo "setting up parameter file"
echo ""
pset hrc_process_events badfile="${OUTDIR}/lev1_hrcs_bad_out.fits"
pset hrc_process_events logfile="${OUTDIR}/lev1_hrcs_out.log"
pset hrc_process_events time_offset=0
pset hrc_process_events instrume=HRC-S
pset hrc_process_events degapfile="${INDIR}/hrcsD1999-11-08gapN0002.fits"
pset hrc_process_events do_ratio=yes
pset hrc_process_events do_amp_sf_cor=yes

pset hrc_process_events badpixfile="${INDIR}/hrcf01246_000N001_bpix1.fits"
pset hrc_process_events acaofffile="${INDIR}/hrcf01246_000N001_aoff1.fits"
pset hrc_process_events alignmentfile="${INDIR}/hrcf01246_000N001_soff1.fits"
pset hrc_process_events gainfile=NONE
pset hrc_process_events ADCfile=NONE
pset hrc_process_events ampsfcorfile="${INDIR}/hrcsAMPSFAug21.fits"
pset hrc_process_events tapfile="NONE"
pset hrc_process_events hypfile="${INDIR}/hrcsD1999-07-23fptestN0001.fits"
pset hrc_process_events ampsatfile="${INDIR}/hrcsD1999-07-23sattestN0001.fits"
pset hrc_process_events evtflatfile="${INDIR}/hrcsD1999-07-22eftestN0001.fits"

pset hrc_process_events cfu1=1.18
pset hrc_process_events cfu2=-0.16
pset hrc_process_events cfv1=1.11
pset hrc_process_events cfv2=-0.1
pset hrc_process_events amp_gain=75
}  # subroutine :  pset_hrcS_low_rangelev
######################################################################
# subroutine
# pset_S_172()
# pre-set parameters for hrc_process_events

pset_S_172()
{
echo "setting up parameter file"
echo ""
pset  hrc_process_events badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits
pset  hrc_process_events  acaofffile=${INDIR}/pcadf052944519N002_asol1.fits
pset  hrc_process_events  alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits
pset  hrc_process_events  gainfile=NONE
pset  hrc_process_events  ADCfile=NONE
pset  hrc_process_events  degapfile="${INDIR}/hrcsD1999-07-22gapN0002.fits"
pset  hrc_process_events  hypfile="${INDIR}/hrcsD1999-07-22fptestN0004.fits"
pset  hrc_process_events ampsfcorfile="${INDIR}/hrcsD1999-07-22amp_sf_corN0001.fits"
pset  hrc_process_events  tapfile="$INDIR/hrcsD1999-07-22tapringN0001.fits"
pset  hrc_process_events  ampsatfile="${INDIR}/hrcsD1999-07-22sattestN0002.fits"
pset  hrc_process_events  evtflatfile="${INDIR}/hrcsD1999-07-22eftestN0001.fits"
pset  hrc_process_events  badfile=${OUTDIR}lev1_bad_evts.qp
pset  hrc_process_events  logfile="${OUTDIR}/lev1_hrcs_out.log"

pset  hrc_process_events  cfu1=1
pset  hrc_process_events  cfu2=0
pset  hrc_process_events  cfv1=1
pset  hrc_process_events  cfv2=0
pset  hrc_process_events  time_offset=0
pset  hrc_process_events  amp_gain=52.9

pset  hrc_process_events  instrume="hrc-s"
pset  hrc_process_events  do_amp_sf_cor=yes
pset  hrc_process_events  do_ratio=yes
} # subroutine : pset_S_172
######################################################################
# Initialization

# !!3
toolname="hrc_process_events"

# set up list of tests
# !!4
alltests="S_warn_nom I_obsfile S_rmNewKey S_addNewKey hrc_I hrc_I_a hrc_I_b hrc_I_c hrc_S hrcS_no_rangelev hrcS_low_rangelev S_no_ampsfcor S_no_ampsfcor2 S_172 S_172_a S_172_a2 S_172_a3 S_172_b S_2582_gainTab S_2582_gainTab_b S_2582_gainTab_c S_1246_gainTab S_1246_gainTab_b S_1246_gainTab_c S_1246_gainTab_d new_I_gainImg new_I_gainImg_b new_I_gainImg_c new_I_gainImg_d"

# "short" test to run
# !!5
shortlist="S_warn_nom I_obsfile S_rmNewKey S_addNewKey hrc_I hrc_S hrcS_no_rangelev hrcS_low_rangelev S_no_ampsfcor S_no_ampsfcor2 S_172 S_172_a"


# compute date string for log file
DT=`date +'%d%b%Y_%T'`


# convenience definitions
OUTDIR=$TESTOUT/$toolname
SAVDIR=$TESTSAV/$toolname
INDIR=$TESTIN/$toolname
LOGDIR=$TESTLOG/$toolname

# set up log file name
LOGFILE=$LOGDIR/${toolname}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${toolname}_log.*



# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi


# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR 
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi

# Make sure we have an input directory
if test -d $INDIR ; then
:
else
  if test $? -ne 0 ; then
    error_exit "can't find input directory : $INDIR"
  fi
fi

# Make sure we have a save directory
if test -d $SAVDIR ; then
:
else
  if test $? -ne 0 ; then
    error_exit "can't find save directory : $SAVDIR"
  fi
fi

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined" 
fi


# check for tools
# if a utility is used in the form "utility <args> > outfile", and 'utility'
# cannot be run, 'outfile' will still be created.  If utility is used on 
# both the output and reference files of a tool the resultant utility output 
# files will both exist and be empty, and will pass a diff.

find_tool dmlist
find_tool dmimgcalc



# announce ourselves
echo ""
echo "${toolname} regression" | tee $LOGFILE
echo "" | tee -a $LOGFILE

# All parameters except verbose should be set anyway, but clear them
# to be safe.
bose=`pget $toolname verbose`
punlearn $toolname
pset $toolname verbose=$bose

# ----- pset geom files -----
# pset geom instruments="$INDIR/telD1999-07-23geomN0002.fits"
# pset geom aimpoints="$INDIR/telD1999-07-23aimptsN0001.fits"
# pset geom tdet="$INDIR/telD1999-07-23tdetN0001.fits"
# pset geom sky="$INDIR/telD1999-07-23skyN0001.fits"


script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do
  # set some parameters; do it inside loop;
  punlearn $toolname
  pset_hpe

  # Init per-test error flag
  # 
  mismatch=1
    
  # delete old outputs
  rm -f $OUTDIR/${testid}*

  # Set up file names
  outfile=$OUTDIR/${testid}.fits
  savfile=$SAVDIR/${testid}.fits

  echo "running $testid" | tee -a $LOGFILE
  echo "" | tee -a $LOGFILE
  ####################################################################
  # run the tool
  case ${testid} in
    #!same as hrc_S, except (obsfile!=NONE) and its NOM keys are diff.
    #!from infile; output data same as hrc_S;
    S_warn_nom)  pset_hrc_S
            test2_string="hrc_process_events \
               infile=${INDIR}/1246_small_evt0a.fits \
               outfile=$outfile \
               obsfile=${INDIR}/1246_obs.par \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test2_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test2_string
            ;;
    #!(11/22/02) same as hrc_S except remove NOM keywords from infile 
    #!NOM won't affect the data if use aoff1; crash test for NOM key; save it;
    S_no_nom)  pset_hrc_S
            pset hrc_process_events infile=${INDIR}/S_no_nom_in.fits
            pset hrc_process_events outfile=$outfile
            pset hrc_process_events badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits
            pset hrc_process_events acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits
            pset hrc_process_events alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits
            test2_string="hrc_process_events \
               infile=${INDIR}/S_no_nom_in.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test2_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test2_string
            ;;
    #!(11/22/02) same as S_172 except remove NOM keywords from infile
    #!NOM affects the data if use asol1; crash test for NOM key; save it;
    S_172_no_nom)  pset_S_172
            pset hrc_process_events infile=${INDIR}/S_172_no_nom_in.fits
            pset hrc_process_events outfile=$outfile
            pset hrc_process_events acaofffile=${INDIR}/pcadf052944519N002_asol1.fits
            pset hrc_process_events alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits
            pset hrc_process_events badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_no_nom_in.fits \
               outfile=$outfile  \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
               acaofffile=${INDIR}/pcadf052944519N002_asol1.fits \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # !! call pset_hrc_I, then Overwrite  obsfile and infile
    # !! test for obsfile!=NONE
    I_obsfile)  pset_hrc_I
            pset hrc_process_events obsfile=$INDIR/axaff00461_001N001_obs0a.par
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a_old.fits \
               outfile=$outfile \
               badpixfile=NONE \
               obsfile=$INDIR/axaff00461_001N001_obs0a.par \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;

    # !! call pset_hrc_S, then Overwrite  obsfile and infile;
    # !! test for obsfile!=NONE
    # !! get keys from obsfile;
    S_rmNewKey)  pset_hrc_S
            test2_string="hrc_process_events \
               infile=${INDIR}/1246_small_evt0a_rmNewKey.fits \
               outfile=$outfile \
               obsfile=${INDIR}/axaff01246_000N001_obs0a.par \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test2_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test2_string
            ;;
    # !! call pset_hrc_S, then Overwrite  obsfile and infile;
    # !! test for obsfile!=NONE
    # !! get keys from evt1 file ;
    S_addNewKey)  pset_hrc_S
            test2_string="hrc_process_events \
               infile=${INDIR}/1246_small_evt0a_addNewKey.fits \
               outfile=$outfile \
               obsfile=${INDIR}/axaff01246_000N001_obs0a.par \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test2_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test2_string
            ;;

    # !!6; gainfile=NONE
    hrc_I)  pset_hrc_I 
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE 
            echo "" | tee -a  $LOGFILE
	    eval  $test1_string
            ;;
    # ref. hrc_I ; set gainfile=hrciD1999-10-04gainN0001.fits
    hrc_I_a)  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=${INDIR}/hrciD1999-10-04gainN0001.fits \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;
    # ref. hrc_I ; set gainfile=CALDB ( hrciD1999-10-04gainN0003.fits )
    hrc_I_b)  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=CALDB \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;
    # ref. hrc_I_b ; set 8 CALDB files ; 
    hrc_I_c)  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               degapfile=CALDB \
               gainfile=CALDB \
               ADCfile=CALDB \
               hypfile=CALDB \
               tapfile=CALDB \
               evtflatfile=CALDB \
               ampsatfile=CALDB \
               ampsfcorfile=CALDB \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;
    # !!7
    hrc_S)  pset_hrc_S
            test2_string="hrc_process_events \
               infile=${INDIR}/1246_small_evt0a.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test2_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test2_string
            ;;

    # !!8
    hrcS_no_rangelev)  pset_hrcS_no_rangelev
            test3_string="hrc_process_events \
               infile=${INDIR}/1246_no_rangelev.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test3_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test3_string
            ;;

    # !!9
    # get the range value from 1246_low_rangelev.par, not from 1246_no_rangelev.fits
    hrcS_low_rangelev)  pset_hrcS_low_rangelev
            test4_string="hrc_process_events \
               infile=${INDIR}/1246_low_rangelev.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
	    echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
	    eval  $test4_string
            ;;

    # !!10
    # get the range value from 1246_low_rangelev.par, not from 1246_no_rangelev.fits
    # use pset_hrcS_low_rangelev to set up the par file, then set do_amp_sf_cor=no;
    S_no_ampsfcor)  pset_hrcS_low_rangelev
            pset hrc_process_events do_amp_sf_cor=no
            test4_string="hrc_process_events \
               infile=${INDIR}/1246_no_rangelev.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    # !! 11
    # get the range value from 1246_small_evt0a.fits, not from 1246_no_rangelev.par
    # use pset_hrcS_low_rangelev to set up the par file, then set do_amp_sf_cor=no;
    S_no_ampsfcor2)  pset_hrcS_low_rangelev
            pset hrc_process_events do_amp_sf_cor=no
            test4_string="hrc_process_events \
               infile=${INDIR}/1246_small_evt0a.fits \
               outfile=$outfile \
               badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
               acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # ( 11/18/02) # obsid_172 : ciao 2.3 hpe issue (ref. points of eqpos[ra,dec])
    # ( 7/2005) test ASOLFILE keyword when acaofffile=@*.lis has two files
    S_172)  pset_S_172
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_in.fits \
               outfile=$outfile  \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
               acaofffile=@${INDIR}/pcadf052944519N002_asol1.lis \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # test degap calib ; test new degap format ; 
    S_172_a)  pset_S_172
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_in.fits \
               outfile=$outfile  \
               degapfile=${INDIR}/hrcsD1999-07-22gaplookupN0001.fits \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
               acaofffile=${INDIR}/pcadf052944519N002_asol1.fits \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # set degapfile=CALDB to test new degap format (hrcsD1999-07-22gaplookupN0002.fits)
    S_172_a2)  pset_S_172
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_in.fits \
               outfile=$outfile  \
               degapfile=CALDB \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
               acaofffile=${INDIR}/pcadf052944519N002_asol1.fits \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # ref: S_172_a2, set 8 CALDB files
    S_172_a3)  pset_S_172
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_in.fits \
               outfile=$outfile  \
               degapfile=CALDB \
               gainfile=CALDB \
               ADCfile=CALDB \
               hypfile=CALDB \
               tapfile=CALDB \
               evtflatfile=CALDB \
               ampsatfile=CALDB \
               ampsfcorfile=CALDB \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
               acaofffile=${INDIR}/pcadf052944519N002_asol1.fits \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # ( 8/2005) test ASOLFILE keyword when acaofffile=/path/a.fits,/path/b.fits
    #           output should be same as S_172.fits
    S_172_b)  pset_S_172
            test4_string="hrc_process_events \
               infile=${INDIR}/S_172_in.fits \
               outfile=$outfile  \
               badpixfile=${INDIR}/hrcf00172_000N002_bpix1.fits \
 acaofffile=${INDIR}/pcadf052944519N002_asol1_A.fits,${INDIR}/pcadf052944519N002_asol1_B.fits \
               alignmentfile=${INDIR}/pcadf052944519N002_asol1.fits > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    #  dph new hrcS new gain table. infile has 9 events only.
    #  no amp_sf corr; no tapring corr;
    S_2582_gainTab ) test4_string="hrc_process_events infile=${INDIR}/S_2582_evt1.fits \
            outfile=$outfile gainfile=${INDIR}/hrcsD1999-08-22tgainN0001.fits \
            degapfile=NONE badpixfile=NONE \
            acaofffile=${INDIR}/S_2582_aoff1.fits \
            alignmentfile=${INDIR}/S_2582_soff1.fits \
            do_ratio=no  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE tapfile=NONE hypfile=NONE \
            ampsatfile=NONE evtflatfile=NONE \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    #  same as S_2582_gainTab except infile has 1000 events.
    #  no amp_sf corr; no tapring corr;
    S_2582_gainTab_b ) test4_string="hrc_process_events infile=${INDIR}/big_S_2582_evt1.fits \
            outfile=$outfile gainfile=${INDIR}/hrcsD1999-08-22tgainN0001.fits \
            degapfile=NONE badpixfile=NONE \
            acaofffile=${INDIR}/S_2582_aoff1.fits \
            alignmentfile=${INDIR}/S_2582_soff1.fits \
            do_ratio=no  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE tapfile=NONE hypfile=NONE \
            ampsatfile=NONE evtflatfile=NONE \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    #  same as S_2582_gainTab_b except use amp_sf corr & use tapring corr; 
    S_2582_gainTab_c ) test4_string="hrc_process_events infile=${INDIR}/big_S_2582_evt1.fits \
            outfile=$outfile gainfile=${INDIR}/hrcsD1999-08-22tgainN0001.fits \
            degapfile=NONE badpixfile=NONE \
            acaofffile=${INDIR}/S_2582_aoff1.fits \
            alignmentfile=${INDIR}/S_2582_soff1.fits \
            do_ratio=no  do_amp_sf_cor=yes \
            ADCfile=NONE ampsfcorfile=${INDIR}/hrcsD1999-07-22amp_sf_corN0001.fits \
            tapfile=$INDIR/hrcsD1999-07-22tapringN0001.fits hypfile=NONE \
            ampsatfile=NONE evtflatfile=NONE \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    #  dph new hrcS gain table. infile is not from archive. 
    #  mjd_obs is added to infile. 
    #  no amp_sf corr ; no tapring corr ;
    S_1246_gainTab ) test4_string="hrc_process_events \
            infile=${INDIR}/upd_S_1246_evt1.fits outfile=${outfile} \
            gainfile=${INDIR}//hrcsD1999-08-22tgainN0001.fits \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE \
            tapfile=NONE hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    # Same as S_1246_gainTab except infile is from archive. 
    # mjd_obs is also added to infile.  Test pi>max status.
    # no amp_sf corr ; no tapring corr ;
    S_1246_gainTab_b ) test4_string="hrc_process_events \
            infile=${INDIR}/arc_S_1246_evt1.fits outfile=${outfile} \
            gainfile=${INDIR}//hrcsD1999-08-22tgainN0001.fits \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE \
            tapfile=NONE hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    # Same as S_1246_gainTab_b except remove mjd_obs from infile(issue WARNING)
    # no amp_sf corr ; no tapring corr ;
    S_1246_gainTab_c ) test4_string="hrc_process_events \
            infile=${INDIR}/arc_S_1246_evt1_N.fits outfile=${outfile} \
            gainfile=${INDIR}//hrcsD1999-08-22tgainN0001.fits \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE \
            tapfile=NONE hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    # Same as S_1246_gainTab_b except use amp_sf corr & use tapring corr ;
    S_1246_gainTab_d ) test4_string="hrc_process_events \
            infile=${INDIR}/arc_S_1246_evt1.fits outfile=${outfile} \
            gainfile=${INDIR}//hrcsD1999-08-22tgainN0001.fits \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=yes \
            ADCfile=NONE ampsfcorfile=${INDIR}/hrcsD1999-07-22amp_sf_corN0001.fits \
            tapfile=$INDIR/hrcsD1999-07-22tapringN0001.fits \
            hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;

    # ref. hrc_I_a; add CHK_GAIN key to gainfile=hrciD1999-10-04gainN0001.fits;
    # use amp_sf corr ; use tapring corr ;
    new_I_gainImg )  pset_hrc_I           
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=${INDIR}/hrciD1999-10-04gain_fake.fits \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;

    # ref. new_I_gainImg; no amp_sf corr ;  no tapring corr ;
    new_I_gainImg_b )  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=${INDIR}/hrciD1999-10-04gain_fake.fits \
               do_amp_sf_cor=no ampsfcorfile=NONE tapfile=NONE \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;

    # ref. new_I_gainImg; use amp_sf corr ; no tapring corr ;
    new_I_gainImg_c )  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=${INDIR}/hrciD1999-10-04gain_fake.fits \
               tapfile=NONE \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;

    # same as new_I_gainImg_b except diff. hrci gain file
    new_I_gainImg_d )  pset_hrc_I
            test1_string="hrc_process_events \
               infile=${INDIR}/3c273_small_evt0a.fits \
               outfile=$outfile badpixfile=NONE \
               gainfile=${INDIR}/hrciD1999-10-04gain_fake2.fits \
               do_amp_sf_cor=no ampsfcorfile=NONE tapfile=NONE \
               acaofffile=${INDIR}/hrcf461_000N001_aoff1.fits \
               alignmentfile=${INDIR}/hrcf461_000N001_soff1.fits > /dev/null"
            echo $test1_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test1_string
            ;;

    # --------------------------- obsolete below ------
    # peter hrcs 3dim gain image. gdropfile is specified.  (obsolete 10/2009)
    2582_gdrop) test4_string="hrc_process_events infile=${INDIR}/S_2582_evt1.fits \
            outfile=$outfile gainfile=${INDIR}/spimeanfits.fits \
            gdropfile=${INDIR}/gaindrop_9.fits \
            degapfile=NONE badpixfile=NONE \
            acaofffile=${INDIR}/S_2582_aoff1.fits \
            alignmentfile=${INDIR}/S_2582_soff1.fits \
            do_ratio=no  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE tapfile=NONE hypfile=NONE \
            ampsatfile=NONE evtflatfile=NONE \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # hrcs new gain file. specify gdropfile but not used. warning.(obsolete 10/2009)
    1246_gainA )
      test4_string="hrc_process_events \
            infile=${INDIR}/S_1246_evt1.fits outfile=${outfile} \
            gainfile=${INDIR}/spimeanfits.fits \
            gdropfile=${INDIR}/gaindrop_1.fits \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE \
            tapfile=NONE hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
    # hrcs new gain file. gdropfile=CALDB (=NONE). warning.(obsolete 10/2009)
    1246_gainB)
      test4_string="hrc_process_events \
            infile=${INDIR}/S_1246_evt1.fits outfile=${outfile} \
            gainfile=${INDIR}/spimeanfits.fits \
            gdropfile=CALDB \
            degapfile=${INDIR}/hrcsD1999-11-08gapN0002.fits \
            badpixfile=${INDIR}/hrcf01246_000N001_bpix1.fits \
            acaofffile=${INDIR}/hrcf01246_000N001_aoff1.fits \
            alignmentfile=${INDIR}/hrcf01246_000N001_soff1.fits \
            do_ratio=yes  do_amp_sf_cor=no \
            ADCfile=NONE ampsfcorfile=NONE \
            tapfile=NONE hypfile=${INDIR}/hrcsD1999-07-23fptestN0001.fits \
            ampsatfile=${INDIR}/hrcsD1999-07-23sattestN0001.fits \
            evtflatfile=${INDIR}/hrcsD1999-07-22eftestN0001.fits \
            cfu1=1.18 cfu2=-0.16 cfv1=1.11 cfv2=-0.1 amp_gain=75 \
            badfile=${OUTDIR}/lev1_hrcs_bad_out.fits \
            logfile=${OUTDIR}/lev1_hrcs_out.log \
            time_offset=0 instrume=HRC-S \
            clob+ verbose=0 > /dev/null"
            echo $test4_string | tee -a  $LOGFILE
            echo "" | tee -a  $LOGFILE
            eval  $test4_string
            ;;
  esac

  ####################################################################
  # if the tool failed to run, set mismatch to 0.
  #
  if test $? -ne 0; then
     echo "" | tee -a $LOGFILE
     echo "$toolname failed to run" | tee -a $LOGFILE
     echo "" | tee -a $LOGFILE
     mismatch=0
  fi


  ####################################################################
  # check the outputs

  # if different tests need different kinds of comparisons, use a 
  #  case ${testid} in...  here

  ####################################################################
  # FITS table    (duplicate for as many tables per test as needed)

  # new output
  # !!10
  #dmlist $outfile header,data,clean > $OUTDIR/${testid}.dmp1  2>>$LOGFILE
  #keyfilter $OUTDIR/${testid}.dmp1 $OUTDIR/${testid}.dmp2  2>>$LOGFILE

  # reference output
  # !!11
  #dmlist $savfile header,data,clean > $OUTDIR/${testid}.dmp1_std  2>>$LOGFILE
  #keyfilter $OUTDIR/${testid}.dmp1_std $OUTDIR/${testid}.dmp2_std \
  #          2>>$LOGFILE

  # compare
  # !!12
  #diff $OUTDIR/${testid}.dmp2 $OUTDIR/${testid}.dmp2_std > \
  #     /dev/null 2>>$LOGFILE
  dmdiff $outfile $savfile tol=$SAVDIR/tolerance > /dev/null 2>>$LOGFILE
  if  test $? -ne 0 ; then
    echo "ERROR: MISMATCH in $outfile" >> $LOGFILE
    mismatch=0
  fi
  ####################################################################
  # FITS image  (duplicate for as many images per test as needed)

  # check image
  # !!13
  # dmimgcalc "$outfile[1]" "$savfile[1]" none tst verbose=0   2>>$LOGFILE
  # if test $? -ne 0; then
  #   echo "ERROR: DATA MISMATCH in $outfile" >> $LOGFILE
  #   mismatch=0
  # fi

  #  Check the header of the image

  # !!14
  # dmlist $outfile header > $OUTDIR/${testid}.dmp1  2>>$LOGFILE
  # keyfilter $OUTDIR/${testid}.dmp1 $OUTDIR/${testid}.dmp2  2>>$LOGFILE

  # !!15
  # dmlist $savfile header > $OUTDIR/${testid}.dmp1_std  2>>$LOGFILE
  # keyfilter $OUTDIR/${testid}.dmp1_std $OUTDIR/${testid}.dmp2_std \
  #            2>>$LOGFILE

  # compare
  # !!16
  # diff $OUTDIR/${testid}.dmp2 $OUTDIR/${testid}.dmp2_std > \
  #       /dev/null 2>>$LOGFILE
  # if  test $? -ne 0 ; then
  #   echo "ERROR: HEADER MISMATCH in $outfile" >> $LOGFILE
  #   mismatch=0
  # fi

  ######################################################################
  # ascii files
  # !!17
  # diff $OUTDIR/${testid}.txt $OUTDIR/${testid}.txt_std > \
  #       /dev/null 2>>$LOGFILE
  # if  test $? -ne 0 ; then
  #   echo "ERROR: TEXT MISMATCH in $OUTDIR/${testid}.txt" >> $LOGFILE
  #   mismatch=0
  # fi

  ####################################################################
  # Did we get an error?
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

  echo "" | tee -a $LOGFILE

done
# end per-test loop
######################################################################


######################################################################
# report results

# blank line
echo "" | tee -a $LOGFILE

if test $script_succeeded -eq 0; then
    echo "${toolname} : PASS" | tee -a $LOGFILE
else
    echo "${toolname} : FAIL" | tee -a $LOGFILE
fi

echo ""
echo "log file in ${LOGFILE}" 
echo ""

exit $script_succeeded
