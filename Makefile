
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

paperM.npy :
	python3 example.py

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
# Full eigen solution
#
paper-eigval.npy : paperM.npy
	python3 eigensolve.py -i paper -o paper

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
# Truncation
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
# Truncation with tapering
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
	--taper-frequency 12.5 \
	--show


#
# Figure 4: Full modal seismogram
#
Figure4:
	python3 plot_seismogram.py -i paperModalSeismogramReplace.npy -c paperSeismogram.npy

#
# Figure 5: Modal seismogram with only eigen values corresponding to less than 12.5Hz
#
Figure5:
	python3 plot_seismogram.py -i paperModalSeismogramReplace_F12.5.npy -c paperSeismogram.npy

#
# Figure 6: Greens functions of full modal solution and truncated from Figures 4 and 5 respectively
Figure6:
	python3 plot_greens.py -i paperModalSeismogramReplace_F12.5.npy -c paperModalSeismogramReplace.npy

#
# Figure 7: Same as Figure 5, but with tapering applied to the Greens function to improve accuracy
#
Figure7:
	python3 plot_seismogram.py -i paperModalSeismogramReplaceTaper_F12.5.npy -c paperSeismogram.npy

#
# Figure A1(a): Modal solution without corrections for Rigid Body Modes
#
FigureA1:
	python3 plot_seismogram.py -i paperModalSeismogram.npy -c paperSeismogram.npy
