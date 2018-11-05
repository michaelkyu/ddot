#!/bin/bash

version="v1.0rc3"

# docker tag ddot-anaconda2:$version michaelkyu/ddot-anaconda2:$version

# docker tag ddot-anaconda3:$version michaelkyu/ddot-anaconda3:$version

docker push michaelkyu/ddot-anaconda2:$version
docker push michaelkyu/ddot-anaconda3:$version
