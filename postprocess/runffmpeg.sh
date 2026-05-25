#!/bin/bash
fmr=5
fmo=5
echo 'Encoding MP4'
ffmpeg -y -framerate $fmr -i rmv%05d.jpg -vf fps=$fmo rmv.mp4
