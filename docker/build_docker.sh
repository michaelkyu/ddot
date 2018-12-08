#!/bin/bash

version="v1.0"

docker build -t michaelkyu/ddot-anaconda2:$version -f Dockerfile-anaconda2 .
docker tag michaelkyu/ddot-anaconda2:$version michaelkyu/ddot-anaconda2:latest

docker build -t michaelkyu/ddot-anaconda3:$version -f Dockerfile-anaconda3 .
docker tag michaelkyu/ddot-anaconda3:$version michaelkyu/ddot-anaconda3:latest

