#!/bin/bash -x

BUILD_PATH=build
NAME=$(basename "${0}" .sh)

LD_LIBRARY_PATH="${BUILD_PATH}:${LD_LIBRARY_PATH}" "${BUILD_PATH}"/test/"${NAME}" \
               > "${BUILD_PATH}"/test/out/"${NAME}".out
diff test/ref/"${NAME}" "${BUILD_PATH}"/test/out/"${NAME}".out
