# cd brainxcan repo
# update 'brainxcan_path' to path to brainxcan repo

module load conda
module load gcc/6.2.0
module load plink/1.90

conda activate /gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/ukb_idp
SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/snakemake/bin/snakemake

$SNMK -s brainxcan/snmk/run.snmk --configfile config.cri_test.yaml -p SBrainXcanAndMR -j5
 
