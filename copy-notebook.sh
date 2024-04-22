#!/bin/bash

if [ ! -f $HOME/python_metadynminer.ipynb ]; then
	cp /usr/local/lib/python_metadynminer.ipynb $HOME
fi
if [ ! -f $HOME/oxytocin ]; then
	cp /usr/local/lib/oxytocin $HOME
fi
if [ ! -f $HOME/tutorial.ipynb ]; then
	cp /usr/local/lib/tutorial.ipynb $HOME
fi
