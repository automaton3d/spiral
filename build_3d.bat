@echo off
rem Build the fused automaton with MSVC/nmake:
rem   spiral_auto.exe      — fused: axis elected from w at the peak
rem The CA core (spiral.c) is compiled headless (-DNO_SDL) and linked
rem with the 3D viewer.
nmake /f Makefile.nmake spiral_auto.exe %*