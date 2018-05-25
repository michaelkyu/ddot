#!/bin/bash

version="v1.0rc0"

docker build -t michaelkyu/ddot-anaconda2:$version -f Dockerfile-anaconda2 .

docker build -t michaelkyu/ddot-anaconda3:$version -f Dockerfile-anaconda3 .

