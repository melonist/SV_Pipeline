#!/usr/bin/env python3

import argparse


class CleanVCF:
    """Call on Lumpy VCF output to remove irregular BND calls, looks for BND called on different chromosomes,
    could be real variants but difficult to interpret"""

    def __init__(self, input_file, output_file):
        self.input = input_file
        self.output = output_file

    def remove_bnd(self):
        out_file = open(self.output, 'w')
        with open(self.input) as vcf:
            line_holder = None
            previous_id = None
            previous_chrom = None
            for line in vcf:
                if '#' in line:
                    out_file.write(line)
                else:
                    line_split = line.split('\t')
                    sv_info = line_split[7].split(';')
                    svtype = sv_info[0].split('=')[1]
                    if svtype != 'BND':
                        out_file.write(line)
                    else:
                        if not line_holder:
                            previous_chrom = line_split[0]
                            previous_id = line_split[2].split('_')[0]
                            line_holder = line
                        else:
                            if previous_chrom == line_split[0] and previous_id == line_split[2].split('_')[0]:
                                out_file.write(line_holder)
                                out_file.write(line)
                                line_holder = None


parser = argparse.ArgumentParser(description='Clean Lumpy .vcf to remove excessive BND')
parser.add_argument('-i', action='store', dest='input',
                    help='Lumpy vcf input, file_name only not complete path')
parser.add_argument('-o', action='store', dest='output',
                    help='Clean lumpy output, file_name only not complete path')
parser.add_argument('-d', action='store', dest='directory',
                    help='path to vcf directory')
args = parser.parse_args()


input_name = args.input
output_name = args.output
directory = args.directory

input_path = directory + input_name
output_path = directory + output_name


fuck_you_bnd = CleanVCF(input_file=input_path, output_file=output_path)
fuck_you_bnd.remove_bnd()
