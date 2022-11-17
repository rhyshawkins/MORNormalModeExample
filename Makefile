
TARGETS = paperM.npy \
	paper-eigval.npy \
	paperSeismogram.npy \
	paperModalSeismogram.npy \
	paperModalSeismogramReplace.npy \
	paperModalSeismogramReplace_F12.5.npy \
	paperModalSeismogramReplaceTaper_F12.5.npy

all : $(TARGETS)

clean :
	rm -f $(TARGETS)


#
# Paper example 5x5 cell Model
#
PAPER_SOURCE_X=0.2e3
PAPER_SOURCE_Y=0.4e3
PAPER_RECEIVER_X=0.8e3
PAPER_RECEIVER_Y=1.0e3
PAPER_DT=2.0e-3
PAPER_NSTEPS=2048
PAPER_T0=0.24
PAPER_F=5.0
PAPER_A=1.0e10

#
# Construct the matrices based on the example mesh and elastic paramters
#
paperM.npy :
	python3 example.py

#
# Compute the seismograms in the time domain. This is the reference solution used to compare
# the accuracy of various modal solutions below
#
paperSeismogram.npy : paperM.npy
	python3 run_lddrk.py -i paper \
	-s $(PAPER_SOURCE_X) -S $(PAPER_SOURCE_Y) \
	-r $(PAPER_RECEIVER_X) -R $(PAPER_RECEIVER_Y) \
	-v \
	-d $(PAPER_DT) \
	-N $(PAPER_NSTEPS) \
	-t $(PAPER_T0) -f $(PAPER_F) -A $(PAPER_A) \
	-o $@

#
# Compute the full eigen solution, required for the various modal seismograms below
#
paper-eigval.npy : paperM.npy
	python3 eigensolve.py -i paper -o paper

#
# Use the full eigen decomposition to compute the seismogram (without correcting rigid-body modes)
#
paperModalSeismogram.npy : paper-eigval.npy paperSeismogram.npy
	python3 run_modal.py -i paper -e paper \
	-s $(PAPER_SOURCE_X) -S $(PAPER_SOURCE_Y) \
	-r $(PAPER_RECEIVER_X) -R $(PAPER_RECEIVER_Y) \
	-v \
	-d $(PAPER_DT) \
	-N $(PAPER_NSTEPS) \
	-t $(PAPER_T0) -f $(PAPER_F) -A $(PAPER_A) \
	-c paperSeismogram.npy \
	-o $@

#
# Same as above, but with rigid body modes corrected
#
paperModalSeismogramReplace.npy : paper-eigval.npy paperSeismogram.npy
	python3 run_modal.py -i paper -e paper \
	-s $(PAPER_SOURCE_X) -S $(PAPER_SOURCE_Y) \
	-r $(PAPER_RECEIVER_X) -R $(PAPER_RECEIVER_Y) \
	-v \
	-d $(PAPER_DT) \
	-N $(PAPER_NSTEPS) \
	-t $(PAPER_T0) -f $(PAPER_F) -A $(PAPER_A) \
	-c paperSeismogram.npy \
	-o $@ \
	--replace-rb

#
# A Truncated modal solution only including frequencies below 12.5 Hz
#
paperModalSeismogramReplace_F12.5.npy : paper-eigval.npy paperSeismogram.npy
	python3 run_modal.py -i paper -e paper \
	-s $(PAPER_SOURCE_X) -S $(PAPER_SOURCE_Y) \
	-r $(PAPER_RECEIVER_X) -R $(PAPER_RECEIVER_Y) \
	-v \
	-d $(PAPER_DT) \
	-N $(PAPER_NSTEPS) \
	-t $(PAPER_T0) -f $(PAPER_F) -A $(PAPER_A) \
	-c paperSeismogram.npy \
	-o $@ \
	--replace-rb \
	-F 12.5 

#
# Same as above, but the Greens' functions are tapered before the convolution step to remove
# ringing artifacts.
#
paperModalSeismogramReplaceTaper_F12.5.npy : paper-eigval.npy paperSeismogram.npy
	python3 run_modal.py -i paper -e paper \
	-s $(PAPER_SOURCE_X) -S $(PAPER_SOURCE_Y) \
	-r $(PAPER_RECEIVER_X) -R $(PAPER_RECEIVER_Y) \
	-v \
	-d $(PAPER_DT) \
	-N $(PAPER_NSTEPS) \
	-t $(PAPER_T0) -f $(PAPER_F) -A $(PAPER_A) \
	-c paperSeismogram.npy \
	-o $@ \
	--replace-rb \
	-F 12.5 -T 80 \
	--taper 0.125 \
	--taper-frequency 12.5 

#
# To replicate figures in the source manuscript
#

#
# Figure 4: Full modal seismogram
#
Figure4: paperModalSeismogramReplace.npy paperSeismogram.npy
	python3 plot_seismogram.py -i paperModalSeismogramReplace.npy -c paperSeismogram.npy

#
# Figure 5: Modal seismogram with only eigen values corresponding to less than 12.5Hz
#
Figure5: paperModalSeismogramReplace_F12.5.npy paperSeismogram.npy
	python3 plot_seismogram.py -i paperModalSeismogramReplace_F12.5.npy -c paperSeismogram.npy

#
# Figure 6: Greens functions of full modal solution and truncated from Figures 4 and 5 respectively
Figure6: paperModalSeismogramReplace_F12.5.npy paperModalSeismogramReplace.npy
	python3 plot_greens.py -i paperModalSeismogramReplace_F12.5.npy -c paperModalSeismogramReplace.npy

#
# Figure 7: Same as Figure 5, but with tapering applied to the Greens function to improve accuracy
#
Figure7: paperModalSeismogramReplaceTaper_F12.5.npy paperSeismogram.npy
	python3 plot_seismogram.py -i paperModalSeismogramReplaceTaper_F12.5.npy -c paperSeismogram.npy

#
# Figure A1(a): Modal solution without corrections for Rigid Body Modes
#
FigureA1: paperModalSeismogram.npy paperSeismogram.npy
	python3 plot_seismogram.py -i paperModalSeismogram.npy -c paperSeismogram.npy
