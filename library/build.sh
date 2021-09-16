for directory in vibron seed Lanczos plot; do
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
    # create lib/
    if [ -d lib ]; then rm -r lib; fi
    mkdir lib
    cd lib
    ln -s ../build/lib$directory.a
    # finish
    cd ../..
done
