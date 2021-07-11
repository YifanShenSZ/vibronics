for directory in v0 v1 v2; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm -f *.exe
    cmake --build .
    cd ../..
done