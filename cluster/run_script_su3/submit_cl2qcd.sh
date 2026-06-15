#!/bin/bash

NODE_CPUS=40

function print_help {
    echo "Usage: ${0} <input_file> <log_file> <jobtime> <partition> <jobname>"
    exit
}

for arg in "${@}"; do
    if [ "${arg}" == '-h' ] || [ "${arg}" == '--help' ]; then
        print_help
    fi
done

if [ $# -ne 5 ]; then
    print_help
fi

# Total memory per node is slightly less than 192 GB
mem_per_cpu=4748
cpus_per_task=${NODE_CPUS}
N_STREAMS=${N_STREAMS:-1}
HOST_SEED_BASE=${HOST_SEED_BASE:-137}
HOST_SEED_STRIDE=${HOST_SEED_STRIDE:-1000003}
GENERATED_INPUT_DIR=${GENERATED_INPUT_DIR:-cluster/run_script_su3/generated}

PROJECT_DIR=${PROJECT_DIR:-/work/mesonqcd/barros/SU23_freezing/Topological-freezing-in-SU2-3-simulations}
jobscript="${PROJECT_DIR}/cluster/run_script_su3/run_cl2qcd_job.sh"

input_file="${1}"
log_file="${2}"

jobtime="${3}"
partition="${4}"
jobname="${5}"
nodes_per_step=1
tasks_per_node=1

submit_one() {
    local submit_input_file="${1}"
    local submit_log_file="${2}"
    local submit_jobname="${3}"

    echo "Submitting job ${submit_jobname} on partition ${partition}, max. time: ${jobtime}."
    echo -e ""
    echo -e "input:\t${submit_input_file}"
    echo -e "log:\t${submit_log_file}"
    echo -e "From input:"
    cat "${submit_input_file}"

    sbatch --partition=${partition} -J"${submit_jobname}" --nodes=${nodes_per_step} --ntasks-per-node=${tasks_per_node} --mem=0 --time=${jobtime} "${jobscript}" "${submit_input_file}" "${submit_log_file}"
}

if [ "${N_STREAMS}" -eq 1 ]; then
    submit_one "${input_file}" "${log_file}" "${jobname}"
    exit
fi

mkdir -p "${GENERATED_INPUT_DIR}"
log_dir=$(dirname "${log_file}")
log_base=$(basename "${log_file}")
log_stem="${log_base%.*}"
log_ext="${log_base##*.}"
if [ "${log_stem}" = "${log_ext}" ]; then
    log_ext="log"
fi
mkdir -p "${log_dir}"

input_base=$(basename "${input_file}")
input_stem="${input_base%.*}"

for stream_id in $(seq 0 $((N_STREAMS - 1))); do
    stream_seed=$((HOST_SEED_BASE + stream_id * HOST_SEED_STRIDE))
    stream_label="s${stream_id}_seed${stream_seed}"
    stream_input_file="${GENERATED_INPUT_DIR}/${input_stem}_${stream_label}.txt"
    stream_log_file="${log_dir}/${log_stem}_${stream_label}.${log_ext}"
    stream_jobname="${jobname}_s${stream_id}"

    awk -F= -v seed="${stream_seed}" -v label="${stream_label}" '
        $1 == "hostSeed" {
            print "hostSeed=" seed
            next
        }
        $1 == "confPrefix" {
            val = $2
            for (i = 3; i <= NF; i++) {
                val = val "=" $i
            }
            sub(/\/configs\/conf\.$/, "_" label "/configs/conf.", val)
            print "confPrefix=" val
            next
        }
        $1 == "PRNGPrefix" {
            val = $2
            for (i = 3; i <= NF; i++) {
                val = val "=" $i
            }
            sub(/\/prng\/prng\.$/, "_" label "/prng/prng.", val)
            print "PRNGPrefix=" val
            next
        }
        $1 == "gaugeObsPrefix" {
            val = $2
            for (i = 3; i <= NF; i++) {
                val = val "=" $i
            }
            sub(/\/gaugeObs$/, "_" label "/gaugeObs", val)
            print "gaugeObsPrefix=" val
            next
        }
        { print }
    ' "${input_file}" > "${stream_input_file}"

    submit_one "${stream_input_file}" "${stream_log_file}" "${stream_jobname}"
done

exit
