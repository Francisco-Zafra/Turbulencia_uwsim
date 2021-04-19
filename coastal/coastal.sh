#!/bin/bash

dist=25

	do

	uwsim -Tgaussian -b10 -d10 -c0.398 -w0.55 -z$dist -a10 -f180 -n1e9 -o "LD_coastal_z"$dist"_hg.rws"

	uwsim -i "LD_coastal_z"$dist"_hg.rws" -a10 -f180 -t0:0.1:3 -o "LD_coastal_z"$dist"_a10_fov180_hg.json"

	done

echo "Finalizado"
