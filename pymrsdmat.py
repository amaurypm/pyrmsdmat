#!/usr/bin/env python3
##
## Superimpose a set of protein structures and report a RSMD matrix, in CSV and Mega-compatible formats, using Pymol as a module.
##
## Amaury Pupo Merino
## amaury.pupo@gmail.com
##
## This script is released under GPL v3.
##

## Importing modules
import argparse
import os
import sys
import numpy as np
from pymol import cmd

## Functions
"""
    rootname

Get root file name, excluding rest of the path and extensions.
"""
def rootname(filename):
    rootname, ext = os.path.splitext(os.path.basename(filename))
    if ext.lower() == ".gz":
        rootname = os.path.splitext(rootname)[0]

    return rootname

def write_csv(filename, unique_files, mat):
	with open(filename, "w") as output_file:
		output_file.write("structures")

		struct_names = [rootname(f) for f in unique_files]

		for struct_name in struct_names:
			output_file.write(",{}".format(struct_name))

		output_file.write("\n")

		for i in range(mat.shape[0]):
			for j in range(mat.shape[1]):
				if j == 0:
					output_file.write("{}".format(struct_names[i]))

				if i > j:
					output_file.write(",{:f}".format(mat[i,j]))

				else:
					output_file.write(",")

			output_file.write("\n")

def write_meg(filename, unique_files, mat):
	with open(filename, "w") as output_file:
		output_file.write("#mega\n")
		output_file.write("!Title: RMSD matrix;\n")
		output_file.write("!Format DataType=Distance DataFormat=LowerLeft NTaxa={:d};\n".format(len(unique_files)))
		output_file.write("!Description\n")
		output_file.write("\tRMSD bewteen structures calculated with pyrmsdmat\n;\n\n")

		struct_names = [rootname(f) for f in unique_files]

		for (i, struct_name) in enumerate(struct_names):
			output_file.write("[{:d}] #{}\n".format(i+1, struct_name))

		output_file.write("\n[     ")

		for i in range(len(struct_names)):
			output_file.write("{:9d}".format(i+1))

		output_file.write("  ]\n")

		for i in range(mat.shape[0]):
			for j in range(mat.shape[1]):
				if j == 0:
					output_file.write("[{:2d}]   ".format(i+1))

				if i > j:
					output_file.write("{:9.3g}".format(mat[i,j]))
				else:
					output_file.write(" " * 9)

			output_file.write("\n")

## Main
def main():
	"""Main function.
	"""
	parser=argparse.ArgumentParser(description="Superimpose a set of protein structures and report a RSMD matrix, in CSV and Mega-compatible formats, using Pymol as a module")
	parser.add_argument('structure', nargs='+', help='Any format supported by Pymol.')
	parser.add_argument('-o', '--output', default='rmsd_matrix', help='output files base name [default: %(default)s].')
	parser.add_argument('-v', '--version', action='version', version='1.0', help="Show program's version number and exit.")

	args=parser.parse_args()

	unique_files = list(set(args.structure))
	unique_files.sort()

	n = len(unique_files)

	if n < 2:
	    parser.error("ERROR: At least two distinct structures are required.")

	rmsd_mat = np.zeros((n, n)) - 1.0 # -1.0 is a nonsense rmsd value.

	for c in range(n-1):
	    try:
	        cmd.load(unique_files[c])

	    except:
	        sys.stderr.write("WARNING: Can not load structure {}. Ignoring it. Corresponding RMSD values will be set to -1.0\n".format(unique_files[c]))
	        continue

	    for r in range(c + 1, n):
	        try:
	            cmd.load(unique_files[r])

	        except:
	            write(sys.stderr, "WARNING: Can not load structure {}. Ignoring it. Corresponding RMSD values will be set to -1.0\n".format(unique_files[c]))
	            continue

	        rmsd_mat[r, c] = cmd.align("{} and name CA".format(rootname(unique_files[r])), "{} and name CA".format(rootname(unique_files[c])))[0]

	        cmd.delete(rootname(unique_files[r]))

	    cmd.delete(rootname(unique_files[c]))

	write_csv(args.output + ".csv", unique_files, rmsd_mat)
	write_meg(args.output + ".meg", unique_files, rmsd_mat)

## Running the script
if __name__ == "__main__":
        main()
