#!/usr/bin/env bash

# Script for create and run a active container of the image build "kinetic_gas:latest" created from "build_container.sh" script.

# Detailed Flags Explanation
# -it : Allocate a terminal inside the container and keep it interactive.
# --rm: It will delete itself if you exit the container.
# --name: name of the active container.
# -v: mount the volume from the host git repo to the /root/code directory inside the active container. This makes it
# possible to edit code on the host and then build the code in the container.
# --privileged: Docker will enable almost all capabilities inside the container as the host.
# $image_name: The name of the image build you utilize to create a container.

image_name="kinetic_gas:latest"
docker run -it --rm \
 --name kinetic_gas_container \
 --privileged \
 -v $(pwd)/:/root/code \
 $image_name



