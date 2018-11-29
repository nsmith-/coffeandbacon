#!/bin/bash

while getopts "b" opt; do
  case $opt in
    b) REBUILD=1;;
  esac
done

if [[ $REBUILD || ! $(docker image ls -q coffea) ]]; then
  pushd container
  docker build --no-cache -t coffea .
  popd
fi

docker run -p 8888:8888 -v $PWD/notebooks:/opt/notebooks coffea
