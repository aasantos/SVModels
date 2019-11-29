#/usr/bin
g++ -c func.cpp -O2 -m64 -o func.o
g++ -c DLM.cpp -O2 -m64 -o DLM.o
#
g++ -c SV.cpp -O2 -m64 -o SV.o
g++ -c SV2F.cpp -O2 -m64 -o SV2F.o
g++ -c SVJ.cpp -m64 -o SVJ.o
g++ -c SVJ2F.cpp -O2 -m64 -o SVJ2F.o
##
g++ -c SVDUR.cpp -O2 -m64 -o SVDUR.o
g++ -c SV2FDUR.cpp -O2 -m64 -o SV2FDUR.o
g++ -c SVJDUR.cpp -O2 -m64 -o SVJDUR.o
g++ -c SVJ2FDUR.cpp -O2 -m64 -o SVJ2FDUR.o
##
g++ -c GammaDuration.cpp -O2 -m64 -o GammaDuration.o
#
g++ -O2 -m64 func.o DLM.o SV.o mainsv.cpp -o SVEstim
g++ -O2 -m64 func.o DLM.o SV.o SV2F.o mainsv2f.cpp -o SV2FEstim
g++ -O2 -m64 func.o DLM.o SV.o SVJ.o mainsvj.cpp -o SVJEstim
g++ -O2 -m64 func.o DLM.o SV.o SVJ.o SV2F.o SVJ2F.o mainsvj2f.cpp -o SVJ2FEstim
#
g++ -O2 -m64 func.o DLM.o SV.o GammaDuration.o SVDUR.o mainsvdur.cpp -o SVDUREstim
g++ -O2 -m64 func.o DLM.o SV.o GammaDuration.o SVDUR.o SVJ.o SVJDUR.o mainsvjdur.cpp -o SVJDUREstim
g++ -O2 -m64 func.o DLM.o SV.o GammaDuration.o SVDUR.o SV2F.o SV2FDUR.o mainsv2fdur.cpp -o SV2FDUREstim
g++ -O2 -m64 func.o DLM.o SV.o SV2F.o SVJ.o SVJ2F.o SVDUR.o SV2FDUR.o SVJDUR.o SVJ2FDUR.o GammaDuration.o mainsvj2fdur.cpp -o SVJ2FDUREstim
#
#
rm func.o
rm DLM.o
#
rm SV.o
rm SV2F.o
rm SVJ.o
rm SVJ2F.o
#
rm SVDUR.o
rm SV2FDUR.o
rm SVJDUR.o
rm SVJ2FDUR.o
#
rm GammaDuration.o 
