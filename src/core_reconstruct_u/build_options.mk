PWD=$(shell pwd)
EXE_NAME=reconstruct_u
NAMELIST_SUFFIX=reconstruct_u
override CPPFLAGS += -DCORE_RECONSTRUCT_U

report_builds:
	@echo "RECONSTRUCT_U=reconstruct_u"
