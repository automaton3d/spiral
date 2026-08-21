@echo off
rem Build the interactive 3D viewer (spiral_3d.exe) with MSVC/nmake.
rem The CA core (spiral.c) is compiled headless (-DNO_SDL) and linked
rem with the 3D viewer (spiral_3d.c).
nmake /f Makefile.nmake spiral_3d.exe %*