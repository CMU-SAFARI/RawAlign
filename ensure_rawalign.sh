#run this script to ensure that the rawalign executable is in PATH
#run as source ensure_rawalign.sh
#this will ensure that rawalign is installed and in PATH
#running bash ensure_rawalign.sh will not work

#check that the rawalign executable exists at bin/rawalign
if [ ! -f bin/rawalign ]; then
    echo "rawalign executable not found at bin/rawalign"
    echo "compiling..."
    make
fi

#ensure that the rawalign executable is in PATH
if [[ ":$PATH:" != *":$PWD/bin:"* ]]; then
    echo "Adding bin to PATH"
    export PATH=$PWD/bin:$PATH
fi
