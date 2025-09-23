#!/bin/bash
#
#SBATCH --job-name=transfer_merged_data
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --partition work1
#SBATCH --output=transfer_merged_data_ex_%j.out
#SBATCH --error=transfer_merged_data_ex_%j.err

#---------------------------------

echo "started working..."

### path/to/the/data
path="/home/aakinad/Phd_data/merged_universal_data.csv"

# Check if the file exists
if [[ ! -f $path ]]; then
    echo "Error: File $path does not exist."
    exit 1
fi


##### start the transfer of data
scp  $path aakinad@secretariat-master.clemson.edu:/data2/gopalan_lab/aakinad/Data

if [[ $? -ne 0 ]]; then
    echo "Error: Data transfer failed."
    exit 1
fi


echo "sending successful..."





