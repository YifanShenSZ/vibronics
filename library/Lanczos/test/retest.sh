echo "Entre symmetry"
cd symmetry
bash retest.sh
cd ..

for directory in Hamiltonian mv Lanczos; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm test.exe
    cmake --build .
    cd ..
    if [ -d input ]; then
        cd input
        ../build/test.exe
        cd ../..
    else
       ./build/test.exe
       cd ..
    fi
done
