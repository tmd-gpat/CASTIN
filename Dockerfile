FROM r-base:4.2.0

RUN apt-get update
RUN apt-get install -y default-jre default-jdk ant vim

RUN ["Rscript", "-e", "install.packages('rJava', repos='https://cloud.r-project.org')"]

ENV R_HOME="/usr/lib/R/"
ENV JRI_DIR="/usr/local/lib/R/site-library/rJava/jri/"

COPY . /repo
WORKDIR /repo

RUN ["ln", "-s", "/usr/local/lib/R/site-library/rJava/jri/JRI.jar", "lib/JRI.jar"]
RUN ant

CMD ["/repo/scripts/run_paired.sh"]
