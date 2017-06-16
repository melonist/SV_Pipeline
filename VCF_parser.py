#! usr/env/ python3

import numpy as np


class MendelianError:

    def __init__(self, vcf_iterator):
        self.vcf = vcf_iterator
        self.mendelian_stats = None

    def mendelian_error_rate(self, mother, father, offspring):
        mendelian_correct = []
        mendelian_error = []
        mother_index = self.vcf.header.index(mother)
        father_index = self.vcf.header.index(father)
        offspring_index = self.vcf.header.index(offspring)
        count = 0
        for line in self.vcf:
            line_strip = line.strip('\n')
            line_split = line_strip.split('\t')
            mother_sample = line_split[mother_index].split(':')
            try:
                mother_genotype = [int(i) for i in mother_sample[0].split('/')]
                father_sample = line_split[father_index].split(':')
                father_genotype = [int(i) for i in father_sample[0].split('/')]
                possible_genotypes = [sorted([mother_genotype[0], father_genotype[0]]),
                                      sorted([mother_genotype[0], father_genotype[1]]),
                                      sorted([mother_genotype[1], father_genotype[0]]),
                                      sorted([mother_genotype[1], father_genotype[1]])]
                offspring_sample = line_split[offspring_index].split(':')
                offspring_genotype = [int(i) for i in offspring_sample[0].split('/')]
                count += 1
                if mother_genotype == '0/0' and father_genotype == '0/0' and offspring_genotype == '0/0':
                    count -= 1
                    continue
                elif offspring_genotype in possible_genotypes:
                    mendelian_correct.append(offspring_genotype)
                else:
                    mendelian_error.append(offspring_genotype)
            except ValueError:
                continue
        self.mendelian_stats = str('Mendelian Error:' + str(round((len(mendelian_error) / count), 3)) +
                                   '\nVariants:' + str(count) + '\nErrors:' + str(len(mendelian_error)))


class vcfIterator(object):

    def __init__(self, file):
        self.file = file
        with open(self.file) as vcf:
            info = []
            alt = []
            format = []
            sample = []
            header = None
            self.count = 0
            for line in vcf:
                if '##INFO' in line:
                    line_strip = line.strip('\n')
                    equal_split = line_strip.split('=')
                    comma_split = equal_split[2].split(',')
                    info.append(comma_split[0])
                elif '##ALT' in line:
                    line_strip = line.strip('\n')
                    equal_split = line_strip.split('=')
                    comma_split = equal_split[2].split(',')
                    alt.append(comma_split[0])
                elif '##FORMAT' in line:
                    line_strip = line.strip('\n')
                    equal_split = line_strip.split('=')
                    comma_split = equal_split[2].split(',')
                    format.append(comma_split[0])
                elif '##SAMPLE' in line:
                    line_strip = line.strip('\n')
                    equal_split = line_strip.split('=')
                    comma_split = equal_split[2].split(',')
                    sample.append(comma_split[0])
                elif '#CHROM' in line:
                    line_strip = line.strip('\n')
                    header = line_strip.split('\t')
                    break
                self.count += 1
            self.info = info
            self.alt = alt
            self.format = format
            self.sample = sample
            self.header = header

    def __iter__(self):
        count = 0
        with open(self.file) as vcf:
            for line in vcf:
                if count > self.count:
                    yield (line)
                count += 1


class VariantStats:

    def __init__(self, vcf_iterator):
        self.vcf = vcf_iterator
        self.stats = None
        self.stats_header = None

    def collect_variants(self):
        # Sort SV by type taken from alt header in lumpy vcf, bin by type, calculate allele frequency, save position
        collection = [[] for i in self.vcf.alt]
        collection.append([])
        for line in self.vcf:
            line_split = line.split('\t')
            sv_info = line_split[self.vcf.header.index('INFO')].split(';')
            svtype = sv_info[0].split('=')[1]
            try:
                svlen = int(sv_info[1].split('=')[1])
            except ValueError:
                svlen = np.NaN
            try:
                position = line_split[0] + ':' + line_split[1] + '-' + sv_info[2].split('=')[1]
                alt_allele_count = 0
                allele_count = 0
            except IndexError:
                postion = 'NA'
                alt_allele_count = 0
                allele_count = 0
            for sample in line_split[9:]:
                sample_gt = sample.split(':')[0]
                try:
                    x = [int(i) for i in sample_gt.split('/')]
                    if x[0] != 0:
                        alt_allele_count += 1
                    if x[1] != 0:
                        alt_allele_count += 1
                    allele_count += 2
                except ValueError:
                    continue
            try:
                allele_frequency = round(alt_allele_count/allele_count, 7)
            except ZeroDivisionError:
                allele_frequency = np.NaN
            variants_stats = [position, svlen, allele_frequency, allele_count]
            try:
                sv_index = self.vcf.alt.index(svtype)
                collection[sv_index].append(variants_stats)
            except ValueError:
                collection[-1].append(variants_stats)
        self.stats = collection
        self.stats_header = self.vcf.alt
