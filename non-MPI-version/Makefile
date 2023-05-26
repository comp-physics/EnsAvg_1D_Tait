ALLFILES = m_globalvar.f90 m_misc.f90 m_inoutvar.f90 m_bubbles.f90 m_rhsvar.f90 m_timesplit.f90 m_rhs.f90 m_timemarch.f90 p_main.f90
COMPILE =gfortran -O2 #-fdefault-real-16
# COMPILE =g95 -O2

all: $(ALLFILES)
	$(COMPILE)  $(ALLFILES) -lm -o b1d

clean: 
	/bin/rm -f *.o *.mod
