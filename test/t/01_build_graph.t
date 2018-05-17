#!/usr/bin/env bash

BASH_TAP_ROOT=../deps/bash-tap
. ../deps/bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for vg

plan tests 1

is $(./build_graph | wc -l) 1 "graph building with the API"

rm -rf build_graph.d build_graph.dSYM