#!/bin/bash

source config.sh

docker_build_working_dir="./docker"
my_dockerfile="$docker_build_working_dir/Dockerfile"
image_name="bmagnuso/crispr:0.0.10"

DOCKER_BUILDKIT=1 \time \
	-v docker build \
	--progress=plain \
	--rm=true \
	-t $image_name \
	--file $my_dockerfile \
	$docker_build_working_dir \
	2>&1 \
| tee build.log
