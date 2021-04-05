for directory in vibron seed Lanczos; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm -f lib*
    cmake --build .
    cd ../..
done