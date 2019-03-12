FROM rocker/tidyverse:3.5.2

RUN apt-get update && apt-get install -y python python-pip virtualenv curl less git zlib1g-dev

RUN Rscript -e "install.packages(c('devtools', 'dplyr'), repos='https://cran.rstudio.com/')" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN Rscript -e 'devtools::install_github("rstudio/shinydashboard")'
RUN Rscript -e 'BiocManager::install("Homo.sapiens", version = "3.8", ask=F)'
RUN Rscript -e 'BiocManager::install("ggbio", version = "3.8", ask=F)'

ADD . /home/rstudio/code
ADD .Rprofile /home/rstudio/.Rprofile
