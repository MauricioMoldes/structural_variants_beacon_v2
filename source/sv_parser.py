#!/usr/bin/env python

"""sv_parse.py: Parses the annotsv annotation of a vcf file and generates a beacon compliant json array file ."""

__author__ = "Arnau Soler, Mauricio Moldes"
__version__ = "0.1"
__maintainer__ = "Arnau Soler, Mauricio Moldes"
__email__ = "arnau.soler@crg.eu, mauricio.moldes@crg.eu"
__status__ = "test"

import csv
import json
from pysam import VariantFile
import logging
import sys
import os

logger = logging.getLogger('incomplete_submissions')

global id
global vcf2bff_hostname
global vcf2bff_filein
global vcf2bff_user
global vcf2bff_ncpuhost
global vcf2bff_fileout
global vcf2bff_cwd
global vcf2bff_projectDir
global info_datasetId
global info_genome
global vcf2bff_version
global gene_dds
global annotation_impact
global aminoacid_changes
global molecular_effects_label
global molecular_effects_id
global variant_quality_filter
global variant_quality_qual
global position_start
global position_end
global position_end_integer
global position_start_integer
global position_refseqId
global position_assemblyId
global case_level_data_zygosity
global case_level_data_zigosity_label
global case_level_data_zigosity_id
global case_level_data_biosampleId
global variation
global variation_variantType
global variation_alternateBases
global variation_referenceBases
global variation_location
global variation_location_interval
global variation_location_interval_start
global variation_location_interval_start_value
global variation_location_interval_start_type
global variation_location_interval_end
global variation_location_interval_end_value
global variation_location_interval_end_type
global variation_location_interval_type
global variation_location_type
global variation_location_sequence_id
global variantInternalId
global identifiers
global identifiers_genomicHGVSId

""" creates variant internal id """


def create_internal_id():
    return "internal_id"


def create_file_in():
    return "vcf2bff_filein"


""" wrapper for generating internal information """


def internal_information():
    create_internal_id()
    create_file_in()


""" Receives Genotype, returns zygosity"""


def case_level_data(genotype):
    zygosity = None

    if genotype == '0/1' or genotype == '0|1' or genotype == '1/0' or genotype == '1|0' or genotype == '0':
        zygosity = "GENO:GENO_0000458"
    if genotype == '1/1' or genotype == '1|1' or genotype == '1':
        zygosity = "GENO:GENO_0000136"

    return zygosity


""" Verifies if multiple samples are present, assign each genotype to sample"""


def multi_sample_case_level_data(case_level_data_biosampleId, sample_genotypes):
    samples = []
    genotypes = []
    if "," in case_level_data_biosampleId:  # uses sample field to verify if multi sample
        samples = case_level_data_biosampleId.split(sep=",")  # split sample

        #for samples, sample_genotypes in:
            #return None  # assign genotype and zygosity per sample


""" reads vcf """


def read_vcf():
    bcf_in = VariantFile("../data/delly.vcf.gz")  # auto-detect input format
    bcf_out = VariantFile('-', 'w', header=bcf_in.header)

    for rec in bcf_in.fetch('1', 10585, 14446418):
        bcf_out.write(rec)


def pre_process_tsv():
    os.system("cat ../data/lumpy.annotated.tsv | cut -f1-15 | sort | uniq > ../data/lumpy.annotated.uniq.tsv")


def read_tsv(data):
    with open('../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.uniq.tsv') as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            position_start_integer = row[2]
            position_end_integer = row[3]
            position_start = row[2]
            position_end = row[3]
            position_refseqId = row[1]
            position_assemblyId = "hg19"
            case_level_data_zigosity_label_raw = row[14]
            case_level_data_zigosity_label = case_level_data_zigosity_label_raw.split(sep=":")[0]
            list_genotypes = case_level_data_zigosity_label_raw.split(sep=":")
            case_level_data_zigosity_id = case_level_data(case_level_data_zigosity_label)
            case_level_data_biosampleId = row[6]

            multi_sample_case_level_data(case_level_data_biosampleId, list_genotypes)

            variation_variant_type = row[5]
            variation_alternate_bases = row[9]
            variation_reference_bases = row[8]
            variation_location_interval_start_value = row[2]
            variation_location_interval_start_type = "Number"
            variation_location_interval_end_value = row[3]
            variation_location_interval_end_type = "Number"
            variation_location_type = "SequenceLocation"
            variation_location_sequence_id = row[5]
            variant_internal_id = row[0]
            identifiers_genomic_hgvs_id = row[5]
            variant_quality_filter = row[11]
            variant_quality_qual = row[10]
            variation_location_interval_type = "TODO"

            data['variantInternalId'] = variant_internal_id
            data['identifiers']['genomicHGVSId'] = identifiers_genomic_hgvs_id
            data['_position']['assemblyId'] = position_assemblyId
            data['_position']['end'] = position_end
            data['_position']['endInteger'] = position_end_integer
            data['_position']['refseqId'] = position_refseqId
            data['_position']['startInteger'] = position_start_integer
            data['_position']['start'] = position_start

            data['variantQuality']['QUAL'] = variant_quality_qual
            data['variantQuality']['FILTER'] = variant_quality_filter

            data['variation']['location']['sequence_id'] = variation_location_sequence_id
            data['variation']['location']['type'] = variation_location_type
            data['variation']['location']['interval']['start']['type'] = variation_location_interval_start_type
            data['variation']['location']['interval']['start']['value'] = variation_location_interval_start_value
            data['variation']['location']['interval']['end']['type'] = variation_location_interval_end_type
            data['variation']['location']['interval']['end']['value'] = variation_location_interval_end_value
            data['variation']['location']['interval']['type'] = variation_location_interval_type

            data['variation']['alternateBases'] = variation_alternate_bases
            data['variation']['referenceBases'] = variation_reference_bases
            data['variation']['variantType'] = variation_variant_type

            # TODO
            # data['caseLevelData']['zygosity']['label'] = case_level_data_zigosity_label
            # data['caseLevelData']['zygosity']['id'] = case_level_data_zigosity_id
            # data['caseLevelData']['biosampleId'] = case_level_data_biosampleId

            write_json(data)


def bff_post_processing():
    os.system("sed -i 's/}{/},{/' ../results/g_variants_sv.json")


""" write bff output in json array format """


def read_json():
    ## Transform from bff to python dictionary
    with open('../data/first_json_template.json') as json_file:
        data = json.load(json_file)  # Load json to python

    return data


""" write bff output in json array format """


def write_json(data):
    # Serializing json
    json_object = json.dumps(data, indent=4)

    # Writing to sample.json
    with open("../results/g_variants_sv.json", "a") as outfile:
        outfile.write(json_object)


if __name__ == "__main__":
    try:
        # configure logging
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s [in %(pathname)s:%(lineno)d]'
        logging.basicConfig(format=log_format)
        # execute main function
        data = read_json()  # READS BFF TEMPLATE
        # pre_process_tsv()
        read_tsv(data)  # POPULATE STUFF
        bff_post_processing()

    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)
