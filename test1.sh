#!/bin/bash

for num in 20000000 100000000 200000000
	do
	for seed in 6764 445455 876533
	do
		for numbuckets in 100000 75000 50000 40000 30000 20000 10000 5000 2500 1000
		do
			sequential_bucketsort $num  $numbuckets $seed
		done
		echo
	done
done
