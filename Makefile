CC          = g++ #icpc
CLINKER     = g++ #icpc

#CFLAGS      =   -Wall -O3 -march=pentium
#CFLAGS      =   -Wall -O3 -xHost
#CFLAGS      = -i-fast -lm  
LIBS        = -lm
DEPEND= makedepend

SRC        = HSMD.c ran_uniform.c system.c boxmuller.c readinput.c initialization.c write.c potential.c sample.c gr.c collisiondynamics.c collisioninfo.c pbc.c minimumratio.c compress.c celldetermine.c neighborcell.c makecell.c updatecell.c escapeinfo.c overlapcheck.c collisionupdate.c
OBJS       = HSMD.o ran_uniform.o system.o boxmuller.o readinput.o initialization.o write.o potential.o sample.o gr.o collisiondynamics.o collisioninfo.o pbc.o minimumratio.o compress.o celldetermine.o neighborcell.o makecell.o updatecell.o escapeinfo.o overlapcheck.o collisionupdate.o
EXECS      = HSMD

default: HSMD

all: $(EXECS)

HSMD:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o HSMD $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

HSMD.o: system.h ran_uniform.h 
ran_uniform.o: system.h ran_uniform.h
boxmuller.o: system.h ran_uniform.h
system.o: system.h
readinput.o: system.h ran_uniform.h
initialization.o: system.h ran_uniform.h
write.o: system.h
potential.o: system.h
sample.o: system.h ran_uniform.h
gr.o: system.h
collisioninfo.o: system.h
collisiondynamics.o: system.h
pbc.o: system.h
minimumratio.o: system.h
compress.o: system.h
celldetermine.o: system.h
neighborcell.o: system.h
makecell.o: system.h
updatecell.o: system.h
escapeinfo.o: system.h
overlapcheck.o: system.h
collisionupdate.o: system.h
