#!/bin/bash
dir=$(dirname ${0})
basedir=${dir}/../..
echo ${dir}
pybot -o NONE -l NONE -r NONE -d ${basedir}/build/test/functional/test-reports -x xunit.xml ${dir}/tests.txt
