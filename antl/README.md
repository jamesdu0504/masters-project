To build:

1) Edit the CONFIG file:
 - SOFTWARE_INC and SOFTWARE_LIB should point to the directories where the
inlude files and, respectively, libraries of GMP and NTL are installed
 - update the CFLAGS to whatever works best for the type of machine you
are building this on.  I find it best to copy the ones that GMP selects
automatically which building it.

2) make appl, to build the libarary and the applications in appl

3) doxygen ANTLdoc-config, to build the documentation in doc
