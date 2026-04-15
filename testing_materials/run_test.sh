#!/bin/usr/env bash

set -e

function main {

	for MA in +1 -1 ; do
	
		for MI in -1 -2 -10 +1 ; do
		
			for IN in -1 -10 +1 ; do
		
				echo "matches=${MA}, mismatches=${MI}, indels=${IN}"
		
				python3 main.py \
					--match ${MA} \
					--mismatch ${MI} \
					--indel ${IN} \
					-i testing_materials/example_data/sequences.fasta \
                    -r testing_materials/example_data/reference.fasta \
					-o testing_materials/example_outputs/alignment_matches${MA}_mismatches${MI}_indels${IN}.csv
			
				wait
			
			done
		
		done
	
	done

}

main

exit