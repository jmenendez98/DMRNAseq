FROM continuumio/anaconda-pkg-build:main

RUN conda update -c base conda && \
    conda install -c bioconda -c conda-forge ont-modkit==0.4.1

CMD ["bash"]