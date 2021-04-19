#!/bin/bash

dist=40


	uwsim -Tgaussian -b10 -d10 -c0.151 -w0.245 -z$dist -a10 -f180 -n1e9 -o "LD_clear_z"$dist"_hg.rws"

	uwsim -i "LD_clear_z"$dist"_hg.rws" -a10 -f180 -t0:0.05:1 -o "LD_clear_z"$dist"_a10_fov180_hg.json"


echo "Finalizado"
