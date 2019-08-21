PWD=$(shell pwd)
EXE_NAME=atmosphere_model
NAMELIST_SUFFIX=atmosphere
override CPPFLAGS += -DCORE_ATMOSPHERE

CODE_GEN_OPTS=--cpf-module
export CODE_GEN_OPTS

report_builds:
	@echo "CORE=atmosphere"
