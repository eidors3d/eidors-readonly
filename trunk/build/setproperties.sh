#!/bin/sh

find ../..  \
        -name \*.m \
     -o -name \*.cpp \
     -o -name \*.html \
     -o -name \*.shtml \
    | grep -v .svn \
    | xargs svn propset svn:keywords "Date Author Id"
