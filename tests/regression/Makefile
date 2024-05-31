ROOT_DIR := $(realpath ../..)
include $(ROOT_DIR)/config.mk

all:
	$(MAKE) -C qsim
	$(MAKE) -C qsim_ck

run-simx:
	$(MAKE) -C qsim run-simx
	$(MAKE) -C qsim_ck run-simx

run-rtlsim:
	$(MAKE) -C qsim run-rtlsim
	$(MAKE) -C qsim_ck run-rtlsim

run-opae:
	$(MAKE) -C qsim run-opae
	$(MAKE) -C qsim_ck run-opae

clean:
	$(MAKE) -C qsim clean
	$(MAKE) -C qsim_ck clean

clean-all:
	$(MAKE) -C qsim clean-all
	$(MAKE) -C qsim_ck clean-all
