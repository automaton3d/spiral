@echo off
rem Build the interactive 3D viewers with MSVC/nmake:
rem   spiral_3d.exe        — base 3D viewer
rem   spiral_3d_bcast.exe  — "grid broadcast" edition
rem   momentum.exe         — momentum CA (pulsating sphere, m = 0)
rem   spiral_auto.exe      — fused: axis elected from w at the peak
rem The CA core (spiral.c) is compiled headless (-DNO_SDL) and linked
rem with each 3D viewer.
nmake /f Makefile.nmake spiral_3d.exe spiral_3d_bcast.exe momentum.exe spiral_auto.exe %*
