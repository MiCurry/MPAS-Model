PWD=$(shell pwd)
EXE_NAME=terrain_model
NAMELIST_SUFFIX=terrain
override CPPFLAGS += -DCORE_TERRAIN

report_builds:
	@echo "CORE=terrain"
