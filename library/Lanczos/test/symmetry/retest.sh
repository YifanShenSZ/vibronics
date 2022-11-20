for directory in C1 Cs-irred1 Cs-irred2; do
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
