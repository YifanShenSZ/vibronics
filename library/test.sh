for directory in vibron seed Lanczos; do
    echo
    echo "Entre "$directory"/test"
    cd $directory/test
    bash test.sh
    cd ..
done

for directory in plot; do
    echo
    echo "Entre "$directory"/test"
    cd $directory/test
    # build
    if [ -d build ]; then rm -r build; fi
    mkdir build
    cd build
    cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort ..
    cmake --build .
    cd ..
    # run
    if [ -d input ]; then
        cd input
        ../build/test.exe
        cd ../../..
    else
       ./build/test.exe
       cd ../..
    fi
done
