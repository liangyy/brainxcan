BXCAN=path-to-brainxcan-repo
OUTDIR=path-to-output-directory
OUTPREFIX=prefix-of-output

module load gcc/6.2.0

# set standalone R dependencies
module load R/3.6.4
export R_LIBS=$R_LIB:~/labshare/softwares/rlib/3.6/

# load plink1.9
module load plink/1.90

# load conda 
module load conda

conda activate brainxcan

snakemake -s $BXCAN/brainxcan/snmk/run.snmk --configfile $BXCAN/config.cri_test.yaml -p SBrainXcanAndMR -j5 --config \
  outdir=$OUTDIR \
  brainxcan_path=$BXCAN \    
  prefix=$OUTPREFIX
