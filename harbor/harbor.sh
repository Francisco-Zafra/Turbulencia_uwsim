#!/bin/bash

dist=10

	do

	uwsim -Tgaussian -b10 -d10 -c2.19 -w0.83 -z$dist -a10 -f180 -n1e9 -o "LD_harbor_z10_hg.rws"

	uwsim -i "LD_harbor_z10_hg.rws" -a10 -f180 -t0:0.1:10 -o "LD_harbor_z10_a10_fov180_hg.json"

	done

echo "Finalizado"
