#!/usr/bin/env bash
export NXF_OPTS="-Xms500M -Xmx20G"
tmux new-session -d -s nf1 './nextflow -trace nextflow run main.nf -config nextflow.config'
