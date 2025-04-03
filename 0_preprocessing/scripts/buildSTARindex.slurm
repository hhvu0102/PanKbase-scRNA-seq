#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=starIndex
#SBATCH --mail-type=FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=16G
#SBATCH --time=12:00:00
#SBATCH --account=scjp99
#SBATCH --partition=standard
#SBATCH --mail-user=vthihong@umich.edu
#SBATCH --output=reports/index.%j.out   # STDOUT output file
#SBATCH --error=reports/index.%j.err    # ERR output file

cd /nfs/turbo/umms-scjp-pank/reference/

module load singularity/3.11.5

singularity exec /scratch/scjp_root/scjp99/vthihong/singularity-cache/porchard-default-star-2.7.10a.img STAR --runMode genomeGenerate --genomeDir star_index_oh100 --genomeFastaFiles /nfs/turbo/umms-scjp-pank/reference/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta --sjdbGTFfile /nfs/turbo/umms-scjp-pank/reference/gencode.v39.annotation.gtf --sjdbOverhang 100 --runThreadN 10

#choose overhang = 100 because for single cell data, I have various read length. Author recommended using default overhang https://github.com/alexdobin/STAR/issues/931
