script_dir=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)
plugin_tgz="$script_dir/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz"
final_plugin_dir="$script_dir/extern/hdf5/lib/plugin"
final_plugin_file="$final_plugin_dir/libvbz_hdf_plugin.so"

if [ ! -f "$final_plugin_file" ]; then
    echo "$final_plugin_file not found locally. Downloading from github..."
    pushd "$script_dir"
    wget -O "$plugin_tgz" "https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz"
    tar xf "$plugin_tgz"
    mkdir -p "extern/hdf5/lib/plugin"
    mv --force "ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin/libvbz_hdf_plugin.so" "extern/hdf5/lib/plugin/libvbz_hdf_plugin.so"
    rm -r "ont-vbz-hdf-plugin-1.0.1-Linux"
    rm "$plugin_tgz"
    popd
fi

export HDF5_PLUGIN_PATH="$final_plugin_dir"
