#!/bin/bash




current_dir=$(pwd)
cd ../Data/default/

if [ "$1" == "" ]; then
	cd $(ls -Art | tail -n 1)
else
	cd "$1"
fi

echo $(pwd)
new_dir=$(pwd)
OUTPUT=$(ls network_vertices_???.csv | tail -n 1)
rm *.png
cd "${current_dir}"

python video.py "${OUTPUT}" $1

cd "${new_dir}"

ffmpeg -framerate 5 -i step_%03d.png video.avi
ffmpeg -framerate 5 -i constraints_step_%03d.png video_constraints.avi

cd "${current_dir}"
