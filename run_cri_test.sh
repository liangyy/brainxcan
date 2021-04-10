BXCAN=path-to-brainxcan-repo
OUTDIR=path-to-output-directory
OUTPREFIX=prefix-of-output

module load conda
module load gcc/6.2.0
module load plink/1.90

conda activate /gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/ukb_idp
SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/snakemake/bin/snakemake

$SNMK -s $BXCAN/brainxcan/snmk/run.snmk --configfile $BXCAN/config.cri_test.yaml -p SBrainXcanAndMR -j5 --config \
  outdir=$OUTDIR \
  brainxcan_path=$BXCAN \    
  prefix=$OUTPREFIX
