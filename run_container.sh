#!/bin/bash

if [[ ! $(docker image ls -q coffea) ]]; then
  pushd container
  docker build -t coffea .
  popd
fi

docker run -p 8888:8888 -v $PWD/notebooks:/opt/notebooks coffea
