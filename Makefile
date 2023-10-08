all:rawalign

rawalign:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src
	mv ./src/rawalign ./bin/

subset:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src subset
	mv ./src/rawalign ./bin/

check_dtw:
	@if [ ! -e bin ] ; then mkdir -p ./bin/ ; fi
	+$(MAKE) -C src check_dtw
	mv ./src/check_dtw ./bin/
	./bin/check_dtw

clean:
	rm -rf bin/
	+$(MAKE) clean -C ./src/

.PHONY: check_dtw
	