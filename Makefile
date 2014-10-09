.PHONY: all
all: common_lib libmfree mygsl mfree fem
	@echo Done!

.PHONY: common_lib
common_lib:
	$(MAKE) -e -C common_lib

.PHONY: libmfree	
libmfree:
	$(MAKE) -e -C libmfree
	
.PHONY: mygsl
mygsl:
	$(MAKE) -e -C mygsl
	
mfree: libmfree mygsl
	$(MAKE) -e -C mfree
	cp mfree/mfre mfre

fem: libmfree
	$(MAKE) -e -C fem
	cp fem/fe fe

	
.PHONY: clean
clean: clean_copied_executables
	$(MAKE) -C common_lib clean
	$(MAKE) -C libmfree clean
	$(MAKE) -C mfree clean
	$(MAKE) -C mygsl clean
	$(MAKE) -C fem clean

.PHONY: clean_copied_executables	
clean_copied_executables:
	-rm -f mfre
	-rm -f fe
	
