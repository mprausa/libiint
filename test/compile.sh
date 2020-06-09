#!/bin/bash

g++ -I../libiint/include -I$HOME/.usr/include -I$HOME/.usr/include/flint -L../libiint -L$HOME/.usr/lib test.cpp -o test -liint -lflint -larb -lgmp
