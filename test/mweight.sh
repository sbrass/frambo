#!/bin/bash -x

NAME=$(basename "${0}" .sh)

test/"${NAME}" > build/test/out/"${NAME}".out
python test/hist.py build/test/out/"${NAME}".out
diff test/ref/"${NAME}".hist build/test/out/"${NAME}".hist
