FROM continuumio/miniconda:4.5.11
RUN conda install meta jupyter numba -y --quiet
RUN conda clean -a
RUN pip install --no-cache-dir histbook vega==1.3.0 vegascope meta pandas

WORKDIR /root/
COPY . /root/
RUN git clone http://cdcvs.fnal.gov/projects/nosql-ldrd && cd nosql-ldrd && git apply ../no_c_tools.patch && python setup.py install
EXPOSE 8888
RUN mkdir /opt/notebooks
CMD /opt/conda/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip='0.0.0.0' --port=8888 --no-browser --allow-root
