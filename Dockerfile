FROM jupyter/base-notebook

RUN mkdir /tmp/build
COPY README.md meta.yaml pyproject.toml /tmp/build/
COPY metadynminer/ /tmp/build/metadynminer/

RUN pip3 install /tmp/build

COPY data/oxytocin /home/jovyan
COPY data/oxytocin /usr/local/lib

COPY data/COLVAR_oxytocin /home/jovyan
COPY data/COLVAR_oxytocin /usr/local/lib

COPY tutorial.ipynb /home/jovyan
COPY tutorial.ipynb /usr/local/lib

COPY python_metadynminer.ipynb /home/jovyan
COPY python_metadynminer.ipynb /usr/local/lib

COPY copy-notebook.sh /usr/local/bin/before-notebook.d

RUN mamba install --yes ipywidgets ipympl 

