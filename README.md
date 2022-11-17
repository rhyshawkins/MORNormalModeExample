# MORNormalModeExample

This repository contains some python scripts used to generate some simple examples and demonstrate some basic
concepts in the publication:

R. Hawkins, M. H. Khalid, K. Smetana, J. Trampert, "Model order reduction for seismic waveform modelling: inspiration from normal modes" submitted to Geophysical Journal International 2022

These example codes require python version 3 and numpy/scipy to be installed.

The easiest way to run these from a unix like environment is to use the Makefile provided

make

This will take around 1 minute (primarily for the eigen decomposition)

Various Figures in the manuscript can be recreated through additional make commands, e.g.

make Figure5
make Figure6
make Figure7
make FigureA1

