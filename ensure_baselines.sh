#run as source ensure_baselines.sh
#this will ensure that baselines are installed and in PATH
#running bash ensure_baselines.sh will not work

#check if all baselines are present
if [ ! -d "baselines" ]; then
    echo "baselines directory not found"
    echo "Installing baselines..."
    bash download_baselines.sh
    echo "Baselines installed"
else
    #check if each baseline exists in baselines/bin
    ALL_FOUND=1
    if [ ! -f "baselines/bin/rawhash" ]; then
        echo "rawhash not found"
        ALL_FOUND=0
    fi
    if [ ! -f "baselines/bin/sigmap" ]; then
        echo "sigmap not found"
        ALL_FOUND=0
    fi
    if [ ! -f "baselines/bin/uncalled" ]; then
        #check if the uncalled executable exists in $PATH
        if ! command -v uncalled &> /dev/null; then
            echo "uncalled not found"
            echo "it might not have been installed correctly by pip, or pip's bin directory might not be in PATH"
            echo "ignoring uncalled for now"
            #ALL_FOUND=0
        fi
    fi
    if [ ! -f "baselines/bin/minimap2" ]; then
        echo "minimap2 not found"
        ALL_FOUND=0
    fi

    if [ $ALL_FOUND -eq 0 ]; then
        echo "Installing baselines..."
        #bash download_baselines.sh
        echo "Baselines installed"
    fi
fi

#ensure baselines/bin is in PATH
if [[ ":$PATH:" != *":$PWD/baselines/bin:"* ]]; then
    echo "Adding baselines/bin to PATH"
    export PATH=$PWD/baselines/bin:$PATH
fi
