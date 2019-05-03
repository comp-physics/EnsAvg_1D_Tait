# makefile for 1D bubble flow parallel computations

SHELL=/bin/sh

.SUFFIXES:  .c .o

OBJ1  = m_globalvar.o mpi_setup.o m_misc.o m_inoutvar.o mpi_transfer.o m_bubbles.o m_rhsvar.o m_rhs.o m_timemarch.o m_timesplit.o p_main.o
#flag1 =  -O3 -w -c -C -mismatch -ieee=full
flag1 =  -O3 -w -c -C
flag2 =  -o
#prog = /home/kando/lam/bin/mpif77
#prog = mpif77
prog = mpif90

all            :  $(OBJ1)
		$(prog) $(OBJ1) $(flag2) b1d

m_globalvar.o  :  m_globalvar.f90
		$(prog) $(flag1) m_globalvar.f90

mpi_setup.o    :  mpi_setup.f90
		$(prog) $(flag1) mpi_setup.f90

m_misc.o       : m_misc.f90
		$(prog) $(flag1) m_misc.f90

m_inoutvar.o   :  m_inoutvar.f90
		$(prog) $(flag1) m_inoutvar.f90

mpi_transfer.o :  mpi_transfer.f90
		$(prog) $(flag1) mpi_transfer.f90

m_bubbles.o    :  m_bubbles.f90
		$(prog) $(flag1) m_bubbles.f90

m_rhsvar.o     :  m_rhsvar.f90
		$(prog) $(flag1) m_rhsvar.f90

m_rhs.o        :  m_rhs.f90
		$(prog) $(flag1) m_rhs.f90

m_timemarch.o  :  m_timemarch.f90
		$(prog) $(flag1) m_timemarch.f90

m_timesplit.o  :  m_timesplit.f90
		$(prog) $(flag1) m_timesplit.f90

p_main.o       :  p_main.f90
		$(prog) $(flag1) p_main.f90

.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod
