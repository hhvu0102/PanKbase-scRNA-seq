#!/usr/bin/env bash
#
# This script make lib config files for every individual with scRNA-seq data so we can run Nextflow preprocessing pipeline

#example
#{
#    "libraries": {
#        "9266-VD-1": {
#            "genome": [
#                "hg38"
#            ],
#            "readgroups": {
#                "133155_a": {
#                    "1": "/scratch/scjp_root/scjp0/shared_data/multiome-workshop/fastq/rna/9266-VD-1-GEX_TCGGCTCTAC-CCGATGGTCT_S9_R1_001.fastq.gz",
#                    "2": "/scratch/scjp_root/scjp0/shared_data/multiome-workshop/fastq/rna/9266-VD-1-GEX_TCGGCTCTAC-CCGATGGTCT_S9_R2_001.fastq.gz"
#                }
#            }
#        }
#    }
#}

# For each readgroup, the '1' fastq file corresponds to the sequencing read including the UMI and the nucleus index (and, for 5' GEX, some cDNA); the '2' fastq file refers to the sequencing read representing only cDNA.


helpFunction()
{
	echo ""
	echo "Example usage: bash makeRNAconfig.bash -t \"GSE142465\" -s SRR10751483 -1 SRR10751483_2.fastq.gz -2 SRR10751483_3.fastq.gz -o rna-library-config/GSE142465_library-config.json"
	echo -e "\t-t STR sTudy name. This should be the name of the directory with fastq files."
	echo -e "\t-s STR Sample name.."
       	echo -e "\t-1 STR List of R1 data names."	
	echo -e "\t-2 STR List of R1 data names."
	echo -e "\t-o STR Output which is file name of the lib-config file."
	exit 1 # Exit script after printing help
}


# read command line
while getopts "?:t:s:1:2:o:" opt
do
	case "$opt" in
		t) study="$OPTARG";;
		s) sample="$OPTARG";;
		1) read1="$OPTARG";;
		2) read2="$OPTARG";;
		o) output="$OPTARG";;
		?) helpFunction ;; # Print helpFunction in case parameter is non-existent
	esac
done

echo "		\"$sample\": {"  >> $output
echo "			\"genome\": ["  >> $output
echo "				\"hg38\""  >> $output
echo "			],"  >> $output
echo "			\"readgroups\": {"  >> $output
echo "				\"${sample}\": {"  >> $output
echo "					\"1\": \"/nfs/turbo/umms-scjp-pank/2_IIDP/0_rawData/$study/$read1\", "  >> $output
echo "					\"2\": \"/nfs/turbo/umms-scjp-pank/2_IIDP/0_rawData/$study/$read2\" "  >> $output
echo "				}"  >> $output
echo "			}"  >> $output
echo "		},"  >> $output

