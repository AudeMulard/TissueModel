#!/bin/bash




current_dir=$(pwd)

cd ../Data/default/

cd $(ls -Art | tail -n 1)

echo $(pwd)
new_dir=$(pwd)
OUTPUT=$(ls network_vertices_?.csv | tail -n 1)

cd "${current_dir}"

python video.py "${OUTPUT}"

cd "${new_dir}"
ls
ffmpeg -framerate 5 -i step_%d.png video.avi

cd "${current_dir}"
