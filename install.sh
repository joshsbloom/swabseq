## Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
## change the permissions of the it.
chmod +x Miniconda3-latest-Linux-x86_64.sh
## Run the installation script
./Miniconda3-latest-Linux-x86_64.sh
##
## To activate conda
bash
git clone https://github.com/joshsbloom/swabseq.git
conda env create -f swabseq.yaml
conda activate swabseq
cd swabseq/code/
bash create_tables.sh