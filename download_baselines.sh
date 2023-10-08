### remove baselines directory if it exists
rm -rf baselines

mkdir -p baselines/bin && cd baselines

### Installing SRA Toolkit
wget -qO- https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.1/sratoolkit.3.0.1-ubuntu64.tar.gz | tar xzv; cp -r ./sratoolkit.3.0.1-ubuntu64/bin/* bin/; rm -rf sratoolkit.3.0.1-ubuntu64

##### Step 1 Compiling the tools
### Cloning and compiling RawHash
# git clone --recursive https://github.com/CMU-SAFARI/RawHash.git rawhash && cd rawhash && make && cp ./bin/rawhash ../bin/ && cd ..

# for reproducibility, use the following commit of RawHash
RAWHASH_COMMIT="--branch e9a56fec18007341231af89b32eec9578b1ba622"
git clone ${RAWHASH_COMMIT} --recursive rawhash && cd rawhash && make && cp ./bin/rawhash ../bin/ && cd ..

### Cloning and compiling Sigmap
# git clone --recursive https://github.com/haowenz/sigmap.git sigmap && cd sigmap && make && cp sigmap ../bin/ && cd ..

# for reproducibility, use the following commit of Sigmap
SIGMAP_COMMIT="--branch c9a40483264c9514587a36555b5af48d3f054f6f"
git clone ${SIGMAP_COMMIT} --recursive https://github.com/haowenz/sigmap.git sigmap && cd sigmap && make && cp sigmap ../bin/ && cd ..

### Cloning and compiling UNCALLED
# git clone --recursive https://github.com/skovaka/UNCALLED.git uncalled && cd uncalled/submods/bwa && git pull origin master && cd ../../ && pip3 install git+https://github.com/nanoporetech/read_until_api.git && sed -i 's/read-until==3\.0\.0/read-until/g' setup.py &&pip3 install git+https://github.com/nanoporetech/read_until_api.git && pip3 install . && cd ..

# for reproducibility, use the following commit of UNCALLED and BWA
UNCALLED_COMMIT="--branch 4c0584ab60a811e74664b0b0d241257d39b967ae"
BWA_COMMIT="&& git checkout 139f68fc4c3747813783a488aef2adc86626b01b"
git clone ${UNCALLED_COMMIT} --recursive uncalled && cd uncalled/submods/bwa && git pull origin master ${BWA_COMMIT} && cd ../../ && pip3 install git+

# Downloading and compiling minimap2
wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2; tar -xf minimap2-2.24.tar.bz2; rm minimap2-2.24.tar.bz2; mv minimap2-2.24 minimap2; cd minimap2 && make && cp minimap2 ../bin/ && cd ..
