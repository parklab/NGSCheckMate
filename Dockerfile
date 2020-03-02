FROM continuumio/miniconda:latest

COPY environment.yml .
RUN apt-get update -y && apt-get install -y zlib1g zlib1g-dev gcc 
RUN conda update -n base -c defaults conda && \
    conda env update --name base --file /environment.yml && \
	rm /environment.yml

WORKDIR /checkmate
COPY patterngenerator/ .
COPY ngscheckmate_fastq-source ./ngscheckmate_fastq-source
COPY patterngenerator ./patterngenerator
RUN cd ngscheckmate_fastq-source && alias cc=$(which gcc) && make
RUN cd patterngenerator && make
RUN chmod +x ngscheckmate_fastq-source/ngscheckmate_fastq && cp ngscheckmate_fastq-source/ngscheckmate_fastq /bin/
RUN chmod +x patterngenerator/patternconverter && cp patterngenerator/patternconverter /bin/

COPY tests/ .
COPY ncm_fastq.py ncm_test.py ncm.py /bin/

RUN chmod +x /bin/ncm_fastq.py && chmod +x /bin/ncm.py
