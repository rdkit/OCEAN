#!/bin/bash


ocean_server_URL=${1}
smi_file=${2}

date_str=`date +%Y_%m_%d__%H_%M`

### YOU CAN MODIFY THESE THINGS

datasource="CHEMBL"
output_file="ocean_report_${date_str}.csv"

# using thresholds, t_threshold is the e_value cutoff for targets; c_threshold is the TC threshold for similar compounds
# using numbers, topn_targets is the number of targets to report; topn_compounds the number of similar compounds

# use one of these:
target_parameter="t_threshold=1E-3"
#target_parameter="topn_targets=10"

# and one of these:
compound_parameter="c_threshold=0.3"
#compound_parameter="topn_compounds=10"

###

### YOU SHOULD NOT MODIFY THE FOLLOWING LINES

base_parameter="datasource=${datasource}&${target_parameter}&${compound_parameter}&"

echo "use ocean_server_URL: ${ocean_server_URL}, datasource: ${datasource}, smi_file: ${smi_file}, output_file: ${output_file}"
report_URL="report"

IFS=$'\n'
tasks=($(cat $smi_file))
tasks_count=${#tasks[@]}

get_header="&print_header=1"
IFS=' '
for (( i=0; i<${tasks_count}; i++ ))
do
    task_array=(${tasks[$i]})
    smi=${task_array[0]}
    id=${task_array[1]}

    echo "Collect OCEAN-Report for ${id} -> ${smi}"

    if [ "${i}" -ne "0" ]
    then
	get_header=""
    fi
    query_parameter="&query_id=${id}&query_smiles=${smi}"
    command="wget \"${ocean_server_URL}/${report_URL}?${base_parameter}${query_parameter}${get_header}\" -q -O - >> ${output_file}"
    echo "${command}"
    wget "${ocean_server_URL}/${report_URL}?${base_parameter}${query_parameter}${get_header}" -q -O - >> ${output_file}

done
