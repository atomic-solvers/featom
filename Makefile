##
# featom
#
# @file
# @version 0.1

app/.f90.o:
	$(FC) $(FFLAGS) -c $<
conv: conv.o
	$(FC) $(FFLAGS) -o $@ $^

include Makefile.dep

# end
