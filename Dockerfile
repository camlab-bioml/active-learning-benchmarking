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

RUN git clone --depth 1 https://github.com/tidyverse/ggplot2.git
