#!/bin/bash

FILE_DIR=$1
PATTERN="*.psam"

find $FILE_DIR -name $PATTERN | sed 's/\.[^.]*$//'
