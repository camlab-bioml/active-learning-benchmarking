FROM bioconductor/bioconductor_docker:RELEASE_3_16-R-4.2.3

ADD scripts/install_R_packages.R /tmp/install_R_packages.R
RUN Rscript /tmp/install_R_packages.R

# Install Python 3
RUN apt-get update && apt-get install -y python3

# Install pip for Python 3
RUN apt-get install -y python3-pip

# Install Python packages
RUN pip3 install numpy pandas sklearn joblib

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN Rscript -e "install.packages('pdftools')"
RUN Rscript -e "install.packages('umap')"
RUN Rscript -e "install.packages('DescTools')"
RUN Rscript -e "install.packages('Seurat')"
RUN Rscript -e "remotes::install_github('vh-d/RPortfolioSimilarity')"

RUN git clone https://github.com/Michael-Geuenich/singleCellNet.git
RUN git clone https://github.com/camlab-bioml/leader.git
