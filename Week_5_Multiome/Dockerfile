FROM satijalab/seurat:latest
ENV DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'R.Version()$version.string'
RUN apt update && apt install -y libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libbz2-dev
RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('stuart-lab/signac', ref = 'develop')"
# RUN R -e "install.packages('Signac')"
RUN R -e "install.packages('TFBSTools')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('EnsDb.Mmusculus.v79')"
RUN R -e "BiocManager::install('chromVAR')"
RUN R -e "BiocManager::install('motifmatchr')"
RUN R -e "remotes::install_github('satijalab/azimuth', ref = 'master')"

RUN pip3 install MACS2

# setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies

