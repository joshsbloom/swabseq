FROM centos:7 AS bcl-build

RUN yum update -y \
    && yum install -y wget unzip git python-devel \
    && yum groupinstall -y "Development Tools"


##### get and install bcl2fastq
RUN cd /tmp \
    && wget https://www.dropbox.com/s/idi0xfu0thurk7q/bcl2fastq2-v2-20-0-linux-x86-64.zip \
    && unzip bcl2fastq2-v2-20-0-linux-x86-64.zip \
    && yum install -y bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
RUN rm /tmp/bcl2fastq2*

#ENTRYPOINT ["bcl2fastq"]
#CMD ["--version"]

# Start with rocker/tidyverse base image
FROM rocker/verse:3.6.3
COPY --from=bcl-build /usr/local/bin/bcl2fastq /usr/local/bin/bcl2fastq

# Install extra *nix utils
# x11, mesa, glu1 are so we can install paletteer
RUN apt-get update \
    && apt-get install -y \
    pigz \
    vim \
    git \
    less \
    curl \
    wget \
    parallel \
    python3-pip \
    bzip2 \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-linux/bs?bt_package=latest" -O /usr/local/bin/bs
RUN chmod +x /usr/local/bin/bs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# R Deps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Python Env

RUN pip3 install pandas
RUN pip3 install openpyxl
RUN pip3 install xlrd
RUN pip3 install python-levenshtein
RUN pip3 install scikit-learn
RUN pip3 install numpy
RUN pip3 install seaborn

RUN install2.r --error \
    ggbeeswarm \
    BiocManager \
    stringdist \
    argparser \
    data.table \
    ggpubr \
    viridis \
    reader 

# Install github packages
#RUN installGithub.r \
#    DavisVaughan/furrr

# Install bioconductor packages
RUN R --slave -e "BiocManager::install(c('savR', 'Rqc', 'ShortRead'))"

RUN echo "basedir=/data/Covid/swabseq" >> /root/.bashrc
RUN echo alias swabseq=\'Rscript --vanilla \"$'basedir'/code/countAmplicons.R\" -r \"$'PWD'/\" -t 32 -l 150000000\' >> /root/.bashrc
WORKDIR /home/rstudio 

#RUN apt-get update --fix-missing && \
#    apt-get install -y \
#    bzip2 \
#    ca-certificates \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/*










#WORKDIR /root
#RUN echo "source /home/jbloom/.bash_aliases" >> /root/.bashrc 
#RUN echo "source /home/jbloom/.bashrc" >> /root/.bashrc 
#RUN echo 'alias swabseq=\'Rscript --vanilla "${basedir}code/countAmplicons.R" -r "${PWD}/" -t 32 -l 150000000\'' >> /root/.bashrc
#RUN echo 'p1="$'basedir'code/countAmplicons.R"' >> /root/.bashrc
#RUN echo 'p2="$'PWD'/"' >> /root/.bashrc
# Install miniconda to /miniconda
#ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
#ENV PATH /opt/conda/bin:$PATH



##### install pip
#RUN cd /tmp \
#    && wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm \
#    && yum install -y epel-release-latest-7.noarch.rpm \
#    && yum install -y python-pip
#RUN pip install --upgrade pip



#===============================================================================
# Fix umask so that everything run at rstudio/user level is group-writable
#RUN echo "umask 0002" >> /home/rstudio/.bashrc

#===============================================================================
# Environment Variables
#ENV starcode_version 1.3

#===============================================================================


#RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
#    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
#    rm ~/miniconda.sh && \
#    /opt/conda/bin/conda clean -tipsy && \
#    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
#    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
#    echo "conda activate base" >> ~/.bashrc && \
#    echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/rstudio/.bashrc && \
#    echo "conda activate base" >> /home/rstudio/.bashrc

#WORKDIR /opt/conda
#RUN conda update -y conda \
#    && conda install -y \
#    -c conda-forge \
#    numpy \
#    seaborn \
#    pandas \
#    scikit-learn \
#    python-levenshtein \
#    openpyxl \
#    xlrd

#RUN pip install primer3-py

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# other utils

# starcode
#WORKDIR /usr/local
#RUN set -u; \
#    wget -O starcode.tar.gz "https://github.com/gui11aume/starcode/archive/${starcode_version}.tar.gz" \
#    && tar xzf starcode.tar.gz \
#    && cd "starcode-${starcode_version}" \
#    && make \
#    && cp starcode /usr/local/bin/starcode \
#    && rm ../starcode.tar.gz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# final touches

#RUN echo "alias ll='ls -Alrth'" >> ~/.bashrc
#RUN echo "alias ll='ls -Alrth'" >> /home/rstudio/.bashrc

