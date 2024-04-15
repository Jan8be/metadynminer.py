FROM jupyter/base-notebook

RUN mkdir /tmp/build
COPY README.md meta.yaml setup.py /tmp/build/
COPY metadynminer/ /tmp/build/metadynminer/

RUN pip3 install /tmp/build

COPY data/oxytocin /home/jovyan
COPY data/oxytocin /usr/local/lib

COPY python_metadynminer.ipynb /home/jovyan
COPY python_metadynminer.ipynb /usr/local/lib
COPY copy-notebook.sh /usr/local/bin/before-notebook.d

RUN mamba install --yes ipywidgets ipympl 

