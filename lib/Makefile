include ../common.mk

all: libnpsf.a

LIBSRCDIR=../$(LIBSRCDIRNAME)

libnpsf.a: $(LIBSRCDIR)/libnpsf.a
	ln $^ $@

$(LIBSRCDIR)/libnpsf.a:
	$(MAKE) -C $(LIBSRCDIR)

clean: 
	@rm -f *.a

.PHONY: all clean
