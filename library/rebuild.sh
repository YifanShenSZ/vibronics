for directory in vibron seed Lanczos plot; do
    echo
    echo "Entre "$directory
    cd $directory/build
    rm lib*
    cmake --build .
    cd ../..
done
