#!/bin/bash -x

NAME=$(basename "${0}" .sh)

test/"${NAME}" > build/test/out/"${NAME}".out
diff test/ref/"${NAME}" build/test/out/"${NAME}".out
