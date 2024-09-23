FROM conda/miniconda3:latest

WORKDIR /app
RUN conda install -y -c bioconda htseq
CMD ["bash"]