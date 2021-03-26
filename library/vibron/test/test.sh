for directory in harmonic vibration wfn; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm -f test.exe
    cmake --build .
    ./test.exe
    cd ../..
done
