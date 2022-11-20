for directory in vibron seed Lanczos; do
    echo
    echo "Entre "$directory"/test"
    cd $directory/test
    bash retest.sh
    cd ..
done

for directory in plot; do
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
