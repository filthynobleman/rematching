% Script for compilation of the Mex functions.
% Please, be sure to check all the paths and library names, or the compilation will fail
mex -r2018a ../src/mex/mex_remesh.cpp -I../include -L../install/rmt/lib -L../install/cut/lib -lMexRMT.lib -lRMT.lib -lcut.lib
mex -r2018a ../src/mex/mex_resample.cpp -I../include -L../install/rmt/lib -L../install/cut/lib -lMexRMT.lib -lRMT.lib -lcut.lib