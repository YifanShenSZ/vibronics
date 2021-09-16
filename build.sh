for directory in v0 v1 v2; do
    echo
    echo "Entre "$directory
    cd $directory
    # build
    if [ -d build ]; then rm -r build; fi
    mkdir build
    cd build
    cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..
    cmake --build .
    cd ..
    # link exe
    if [ -f vibronics.exe ]; then rm vibronics.exe; fi
    ln -s build/$directory.exe vibronics.exe
    # finish
    cd ..
done
