#!/bin/sh

awk '/^# analysis/ {filename=$3;next} {print $0 >> filename".tsv"}' $1
