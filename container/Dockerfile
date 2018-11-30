FROM continuumio/miniconda:4.5.11
RUN conda install meta jupyter numba lz4 cloudpickle -y --quiet
RUN conda clean -a
RUN pip install --no-cache-dir histbook vega==1.3.0 vegascope meta pandas uproot

WORKDIR /root/
COPY . /root/
RUN git clone http://cdcvs.fnal.gov/projects/nosql-ldrd && cd nosql-ldrd && git apply ../no_c_tools.patch && python setup.py install && cd ..
# Will be published to pypi.org eventually
RUN git clone https://github.com/CoffeaTeam/fnal-column-analysis-tools.git && cd fnal-column-analysis-tools && python setup.py install && cd ..

EXPOSE 8888
RUN mkdir /opt/notebooks
CMD /opt/conda/bin/jupyter notebook --notebook-dir=/opt/notebooks --ip='0.0.0.0' --port=8888 --no-browser --allow-root
