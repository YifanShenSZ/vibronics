for directory in v0 v1 v2; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm *.exe
    cmake --build .
    cd ../..
done
