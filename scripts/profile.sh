#!/bin/bash

OUTPUT="prof"
TESTS="$(ls test/*.f90 | \
             xargs basename -s .f90 | \
             tr '\n' ' ')"

mkdir -p "${OUTPUT}"

for t in ${TESTS}; do
    "build/test/${t}"
    gprof "build/test/${t}" | \
        gprof2dot -n 0 -e 0 | \
        dot -Tpng -o "${OUTPUT}/${t}.png"
done
