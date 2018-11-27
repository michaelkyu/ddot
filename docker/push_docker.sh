#!/bin/bash

version="v1.0rc3"

docker push michaelkyu/ddot-anaconda2:$version
docker push michaelkyu/ddot-anaconda3:$version

docker push michaelkyu/ddot-anaconda2:$latest
docker push michaelkyu/ddot-anaconda3:$latest
