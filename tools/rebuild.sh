for directory in eig plot2D seed soc; do
    echo
    echo "Entre "$directory
    cd $directory/build
    cmake --build .
    cd ../..
done
