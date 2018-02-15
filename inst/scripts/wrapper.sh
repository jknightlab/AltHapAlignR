#!/bin/bash


# Write here the path to your virtualenv to activate it
if [ ! -z "$1" ]
then
    source "$1/bin/activate"
fi
shift

exec "$@"
