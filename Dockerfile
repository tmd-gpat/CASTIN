FROM r-base:4.2.0

RUN apt-get update
RUN apt-get install -y default-jre default-jdk ant

RUN ["Rscript", "-e", "install.packages('rJava', repos='https://cloud.r-project.org')"]

ENV R_HOME="/usr/local/lib/R/"
ENV JRI_DIR="/usr/local/lib/R/site-library/rJava/jri/F"

COPY . /repo
WORKDIR /repo

RUN ["ln", "-s", "/usr/local/lib/R/site-library/rJava/jri/JRI.jar", "lib/JRI.jar"]
RUN ant

ENTRYPOINT ["/repo/scripts/run_paired.sh"]
