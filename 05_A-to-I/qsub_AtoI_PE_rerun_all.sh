#!/bin/bash
set -euo pipefail
qsub 03_pipeline_G21-4_0515.pbs
qsub 03_pipeline_G21-4_0305.pbs
qsub 03_pipeline_G21-7_0322.pbs
