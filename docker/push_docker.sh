#!/bin/bash

version="v0.2rc1"

docker tag ddot-anaconda2:$version michaelkyu/ddot-anaconda2:$version

docker tag ddot-anaconda3:$version michaelkyu/ddot-anaconda3:$version
