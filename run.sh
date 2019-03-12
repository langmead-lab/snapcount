#!/bin/sh

docker run --name snapr --rm \
    -e PASSWORD=snapr \
    -p 8787:8787 \
    -d \
    snapr
