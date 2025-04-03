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
	echo "Example usage: bash makeRNAconfig.bash -s \"HPAP-020\" -d \"HPAP-020_10xscRNA_73876.1803.L005.R1_fastq-data.fastq.gz HPAP-020_10xscRNA_73876.1837.L005.R1_fastq-data.fastq.gz HPAP-020_10xscRNA_73876.L003.R1_fastq-data.fastq.gz HPAP-020_10xscRNA_73876.L006.R1_fastq-data.fastq.gz HPAP-020_10xscRNA_73876.L007.R1_fastq-data.fastq.gz HPAP-020_10xscRNA_73876.L008.R1_fastq-data.fastq.gz\" -i 6 -o rna-library-config/HPAP-020_library-config.json"
	echo -e "\t-s STR Individual name. This should be the name of the directory with fastq files."
       	echo -e "\t-d STR List of R1 data names, protected in quotation marks."	
	echo -e "\t-i INT Number of read groups of an individual."
	echo -e "\t-o STR Output which is file name of the lib-config file."
	exit 1 # Exit script after printing help
}


# read command line
while getopts "?:s:d:i:o:" opt
do
	case "$opt" in
		s) sample="$OPTARG";;
		d) data="$OPTARG";;
		i) number="$OPTARG";;
		o) output="$OPTARG";;
		?) helpFunction ;; # Print helpFunction in case parameter is non-existent
	esac
done

arr=($data)
echo "		\"$sample\": {"  >> $output
echo "			\"genome\": ["  >> $output
echo "				\"hg38\""  >> $output
echo "			],"  >> $output
echo "			\"readgroups\": {"  >> $output
for (( i=1; i<=$number; i++ )); do
	if [ $i -lt $number ]
	then
		echo "				\"${sample}_${i}\": {"  >> $output
		echo "					\"1\": \"/nfs/turbo/umms-scjp-pank/1_HPAP/data/$sample/data/${arr[$i-1]}\", "  >> $output ### change /nfs/turbo/umms-scjp-pank/1_HPAP/data to directory where data is stored
		r2=`echo ${arr[$i-1]} | sed 's/R1/R2/g'`
		echo "                                  	\"2\": \"/nfs/turbo/umms-scjp-pank/1_HPAP/data/$sample/data/$r2\" "  >> $output
		echo "				},"  >> $output
	else
		echo "                          	\"${sample}_${i}\": {"  >> $output
		echo "                                  	\"1\": \"/nfs/turbo/umms-scjp-pank/1_HPAP/data/$sample/data/${arr[$i-1]}\", "  >> $output
		r2=`echo ${arr[$i-1]} | sed 's/R1/R2/g'`
		echo "                                  	\"2\": \"/nfs/turbo/umms-scjp-pank/1_HPAP/data/$sample/data/$r2\" "  >> $output
		echo "                          	}"  >> $output
	fi
done
echo "			}"  >> $output
echo "		},"  >> $output


