LDLIBS=-lcfitsio

all: tablist

tablist:: tablist.c gftag_normals.c gftag_normals.h

test_%: %.c %.h
	$(LINK.c) $< $(LOADLIBES) $(LDLIBS) -o $@ -D__DO_MAIN_TEST__
