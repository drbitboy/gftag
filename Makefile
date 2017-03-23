CSPICEtop=$(shell ls -d $(HOME)/cspice 2>/dev/null || ls -d /usr/local/spice/cspice || ls -d /usr/local/cspice || echo /CSPICE_UNKNOWN)
CSPICEinc=$(CSPICEtop)/include
CSPICElib=$(CSPICEtop)/lib/cspice.a

CPPFLAGS=-I$(CSPICEinc)
LDLIBS=$(CSPICElib) -lcfitsio -lm

EXES=gftag tablist test_gftag_normals

all: $(EXES)

tablist:: tablist.c

gftag:: gftag.c gftag_normals.c gftag_normals.h

test_%: %.c %.h
	$(LINK.c) -D__DO_MAIN_TEST__ $< $(LOADLIBES) $(LDLIBS) -o $@

clean:
	$(RM) $(EXES)
