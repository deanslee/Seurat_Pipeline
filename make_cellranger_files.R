sample = 'VMPO_SAL' # sample/family name
runid = 'NB502063_0242' # run id from basespace url
n.samples = 1 # number of samples to demultiplex
indices = c('SI-GA-F9')
samples = c('SAL')
setwd('n/regal/dulac_lab/dlee') # where to save locally
flowcells = c('AHV7LNBGX7')
fastq.names = c('VMPO_SAL') # names of folders where fastqs are stored (don't include '_fastq' in these names)
regal.dir = '/n/regal/dulac_lab/dlee'

#######################################

# download raw files
fileConn<-file(paste('download_',sample,'.sh',sep = ''))
writeLines(paste(
'#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -c 1
#SBATCH -t 0-24:00 # Runtime in D-HH:MM # 1 sample took ~6 hrs to download
#SBATCH -p general
#SBATCH --mem 32000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deanlee@fas.harvard.edu
#SBATCH -o ',regal.dir,'/raw/',sample,'.out # Standard out goes to this file
#SBATCH -e ',regal.dir,'/raw/',sample,'.err # Standard err goes to this file
#SBATCH --open-mode=truncate # Append to log files (useful when jobs are killed and restarted)
#SBATCH --constraint=holyib # Make sure nodes are attached to regal so that they are fast

# Load python
source new-modules.sh
module load python/2.7.14-fasrc02

# Change directory
cd ',regal.dir,'/raw

# Download
# make sure the python script is already in the directory ',regal.dir,'/raw
python BaseSpaceRunDownloader_v2.py -r ',runid,' -a 6d0c7c45ef3542a88377665df9ee20ea', sep=''),fileConn)
close(fileConn)
# submit the resulting .sh file to run on the cluster

#######################################

# make fastq file
csv.table = data.frame(Lane = rep(1:4, times = n.samples), Sample= rep(samples, each = 4) , Index = rep(indices, each = 4))
write.csv(csv.table, file = paste(sample,'.csv', sep = ''), row.names = FALSE)
# move this .csv file to the correct directory

fileConn<-file(paste('mkfastq_',sample,'.sh',sep = ''))
writeLines(paste(
'#!/bin/bash
#SBATCH -N 1 # ensures that all cores are on one machine
#SBATCH -n 16
#SBATCH -t 0-24:00 # runtime in D-HH:MM # 1 sample took ~5 hrs to download
#SBATCH -p general
#SBATCH --mem 128000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deanlee@fas.harvard.edu
#SBATCH -o ',regal.dir,'/fastq/',sample,'.out # Standard out goes to this file
#SBATCH -e ',regal.dir,'/fastq/',sample,'.err # Standard err goes to this file
#SBATCH --open-mode=truncate # Append to log files (useful when jobs are killed and restarted)
#SBATCH --constraint=holyib # Make sure nodes are attached to regal so that they are fast

# Load cellranger
source new-modules.sh
module load cellranger/2.1.0-fasrc01

# Change directory
cd ',regal.dir,'/fastq

# Call cellranger
cellranger mkfastq --run=',regal.dir,'/raw/',runid,' --id=',sample,'_fastq --csv=',regal.dir,'/fastq/',sample,'.csv --localmem=128 --localcores=16', sep=''),fileConn)
close(fileConn)
# submit the resulting .sh file to run on the cluster

#######################################

# make counts file
fastq.loc = vector() # make empty vector to store fastq locations
for (i in 1:length(flowcells)){
  for(name in samples){
    a = paste(regal.dir,'/fastq/',fastq.names[i],'_fastq/outs/fastq_path/',flowcells[i],'/',name,sep='')
    fastq.loc[name] = paste(a,fastq.loc[name], sep =',', collapse = ',')
  }
}
fastq.loc = gsub(',NA','',fastq.loc) # remove NA that was introduced on first round of loop

count.call = vector() # make empty vector to store commands needed to run cellranger count
for (i in 1:length(samples)){
  a=paste('cellranger count --id=',samples[i],' --fastqs=',fastq.loc[i],' --sample=',samples[i],' --transcriptome=',regal.dir,'/counts/mm10-1.2.0_premrna --localcores=24 --localmem=192 --nosecondary --chemistry=SC3Pv2', sep='')
  count.call[i]= a
}

fileConn<-file(paste('count_',sample,'.sh',sep = ''))
writeLines(paste(
'#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -c 24
#SBATCH -t 2-00:00 # runtime in D-HH:MM # 1 sample took ~16 hrs to download
#SBATCH -p general # Partition to submit to
#SBATCH --mem 192000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=deanlee@fas.harvard.edu
#SBATCH -o ',regal.dir,'/counts/',sample,'.out # Standard out goes to this file
#SBATCH -e ',regal.dir,'/counts/',sample,'.err # Standard err goes to this file
#SBATCH --open-mode=truncate # Append to log files (useful when jobs are killed and restarted)
#SBATCH --constraint=holyib # Make sure nodes are attached to regal so that they are fast

# Load cellranger
source new-modules.sh
module load cellranger/2.1.0-fasrc01

# Change directory
cd ',regal.dir,'/counts

# Call cellranger
',print(paste(count.call, collapse='\n')), sep=''),fileConn)
close(fileConn)
# submit the resulting .sh file to run on the cluster
