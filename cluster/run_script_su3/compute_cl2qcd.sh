#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: ${0} <input_file> <log_file>"
    exit
fi

program="/home/mesonqcd/public/barros/cl2qcd_obc_barros/build/su3heatbath"

input_file="${1}"
log_file="${2}"
log_dir="$(dirname "${log_file}")"
run_program="${log_dir}/su3heatbath"

echo -e "input:\t${input_file}"
echo -e "log:\t${log_file}"
echo "From input:"
cat "${input_file}"

mkdir -p "${log_dir}"
for output_key in confPrefix PRNGPrefix gaugeObsPrefix; do
    output_prefix=$(awk -F= -v key="${output_key}" '
        $1 == key {
            val = $2
            for (i = 3; i <= NF; i++) {
                val = val "=" $i
            }
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", val)
            print val
            exit
        }
    ' "${input_file}")
    if [ -n "${output_prefix}" ]; then
        mkdir -p "$(dirname "${output_prefix}")"
    fi
done
ln -sfn "${program}" "${run_program}"

"${run_program}" --inputFile "${input_file}" &> "${log_file}"
