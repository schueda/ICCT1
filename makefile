    CC     = gcc -g
    CFLAGS = 
    LFLAGS = -lm

      PROG = cgSolver
      OBJS = sislin.o \
			 metodos.o \
			 utils.o \
             $(PROG).o

.PHONY: limpa faxina clean distclean purge all

%.o: %.c %.h sislin.h metodos.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG):  $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)


purge:
	@rm -f *~ *.bak
	@rm -f *.o core a.out
	@rm -f $(PROG)

clean: purge

