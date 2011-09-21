.PHONY: all clean tools libs
all:   tools libs
clean:
	$(MAKE) -C tools clean
	$(MAKE) -C libs clean
tools:	libs
	$(MAKE) -C tools all
libs:
	$(MAKE) -C libs all
