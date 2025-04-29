#!/bin/bash

# example docker run

container_name="crispr-pipeline"
base_dir=$(pwd)
docker_image="bmagnuso/crispr:0.0.4"
input_dir="input"
output_dir="output"
docker run --rm -itd --name $container_name -v $base_dir:/repo -v $base_dir/$output_dir:/output -v $base_dir/$input_dir:/input -w /repo $docker_image bash
docker attach $container_name
