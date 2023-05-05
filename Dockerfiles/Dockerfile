# Importing minimal installation binary of ubuntu 20
FROM ubuntu:20.04

# Specify shell script. However, might not be necessarry since bash is default.
SHELL ["/bin/bash", "-c"]

# ESSENTIAL
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -yq apt-utils dialog
RUN apt-get install -yq build-essential software-properties-common
RUN apt-get update

# POTENTIAL DEPENDENCIES
RUN apt-get install -yq git curl vim
RUN apt-get install -yq gcc g++ cmake make # C++ dependencies

# Python3.9 and pip3 installation.
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get install -y python3.9 python3.9-tk python3.9-dev
RUN apt-get -y install python3-pip
RUN /usr/bin/python3.9 -m pip install --upgrade pip
RUN pip3 install --upgrade pip
RUN apt-get -y update

# Defining HOME_DIR as environment variable
ENV HOME_DIR /root
RUN mkdir -p $HOME_DIR/code
ENV CODE_DIR $HOME_DIR/code

# Defining the work directory.
WORKDIR $CODE_DIR

# Copy all files from host to container
COPY .. .
RUN echo "source /root/code/common_scripts.sh" >> ~/.bashrc
RUN pip3 install -r requirements.txt

RUN bash cpp/build_kingas.sh
# RUN python3.9 -m pykingas -test




