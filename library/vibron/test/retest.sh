for directory in harmonic vibration wfn; do
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
