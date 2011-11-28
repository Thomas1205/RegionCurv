include Makefile.common

all: $(LIB) $(OPTDIR) $(DEBUGDIR) $(LIB)/commonlib.debug $(LIB)/commonlib.opt $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o

$(LIB)/commonlib.debug: $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(OPTDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o $(DEBUGDIR)/makros.o $(DEBUGDIR)/combinatoric.o $(DEBUGDIR)/storage1D.o $(DEBUGDIR)/storage2D.o
	ar rs $@ $(DEBUGDIR)/fileio.o $(DEBUGDIR)/stringprocessing.o $(DEBUGDIR)/application.o $(DEBUGDIR)/timing.o $(OPTDIR)/vector.o $(DEBUGDIR)/matrix.o $(DEBUGDIR)/tensor.o  $(DEBUGDIR)/makros.o  $(DEBUGDIR)/combinatoric.o $(DEBUGDIR)/storage1D.o $(DEBUGDIR)/storage2D.o

$(LIB)/commonlib.opt: $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o $(OPTDIR)/storage1D.o $(OPTDIR)/storage2D.o
	ar rs $@ $(OPTDIR)/fileio.o $(OPTDIR)/stringprocessing.o $(OPTDIR)/application.o $(OPTDIR)/timing.o $(OPTDIR)/vector.o $(OPTDIR)/matrix.o $(OPTDIR)/tensor.o  $(OPTDIR)/makros.o  $(OPTDIR)/combinatoric.o $(OPTDIR)/storage1D.o $(OPTDIR)/storage2D.o

include Makefile.finish

