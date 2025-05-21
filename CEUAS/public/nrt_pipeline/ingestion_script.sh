#!/bin/bash
mamba activate cads-obs
export CADSOBS_LOGGING_LEVEL="INFO"
export CADSOBS_LOGGING_FORMAT="CONSOLE"
export CLI_DEBUG=true

dataset="insitu-comprehensive-upper-air-observation-network"

nohup cadsobs make-production --dataset insitu-comprehensive-upper-air-observation-network --start-year YYYY --start-month MM --end-year YYYY --end-month MM --source CUON --config /data/private/config/cdsobs/cdsobs_config.yml >& uv_make_prod_up_to_2010.log &
mkdir forms_jsons/insitu-comprehensive-upper-air-observation-network
cadsobs get_forms_jsons --dataset insitu-comprehensive-upper-air-observation-network --output-dir forms_jsons/insitu-comprehensive-upper-air-observation-network --upload
