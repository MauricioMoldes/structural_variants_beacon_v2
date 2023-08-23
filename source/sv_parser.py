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
import socket
import argparse

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

""" creates hostname """


def create_hostname():
    info_vcf2bff_hostname = socket.gethostname()
    return info_vcf2bff_hostname


""" creates variant file in """


def create_filein(filename):
    info_vcf2bff_filein = os.path.basename(filename)
    return info_vcf2bff_filein


""" creates user name """


def create_user():
    info_vcf2bff_user = os.getlogin()
    return info_vcf2bff_user


""" creates host number cpu """


def create_ncpuhost():
    info_vcf2bff_ncpuhost = os.cpu_count()
    return info_vcf2bff_ncpuhost


""" creates bff file out """


def create_fileout():
    info_vcf2bff_fileout = 'genomicVariantsVcf_test.json.gz'  # Needs modification depending on what we want to put
    return info_vcf2bff_fileout


""" creates path (cwd) to variants file """


def create_cwd():
    info_vcf2bff_cwd = '/fullpath/to/data_test/beacon_4524514test/vcf'  # Needs modification depending on what we want to put
    return info_vcf2bff_cwd


""" creates variants project directory name"""


def create_projectDir():
    info_vcf2bff_projectDir = 'beacon_4524514test'  # Needs modification depending on what we want to put
    return info_vcf2bff_projectDir


""" creates version of the beacon"""


def create_version():
    info_vcf2bff_version = '1.0.0'  # Needs modification depending on what we want to put (come from config file from deploy?)
    return info_vcf2bff_version


""" creates datasetId (default one if none is provided?)"""


def create_datasetId():
    info_datasetId = 'default_beacon_1'  # Needs modification depending on what we want to put
    return info_datasetId


""" creates genome"""


def create_genome():
    info_genome = 'hg19'  # Needs to be automatised (from vcf or tsv?)
    return info_genome


""" wrapper for generating internal information """


def internal_information(data):
    filename = '../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.uniq.tsv'  # File name that was processed (vcf or tsv?)
    data['_info']['vcf2bff']['hostname'] = create_hostname()  # Hostname
    data['_info']['vcf2bff']['filein'] = create_filein(filename)  # File in
    data['_info']['vcf2bff']['user'] = create_user()  # User
    data['_info']['vcf2bff']['ncpuhost'] = create_ncpuhost()  # Number of cpu's
    data['_info']['vcf2bff']['fileout'] = create_fileout()  # File out
    data['_info']['vcf2bff']['cwd'] = create_cwd()  # Path to output folder
    data['_info']['vcf2bff']['cwd'] = create_cwd()  # Path to output folder
    data['_info']['vcf2bff']['projectDir'] = create_projectDir()  # Path to output final project directory
    data['_info']['vcf2bff']['version'] = create_version()  # Version of the Beacon


""" Receives Genotype, returns zygosity"""


def case_level_data(genotype):
    zygosity = None

    if genotype == '0/1' or genotype == '0|1' or genotype == '1/0' or genotype == '1|0' or genotype == '0':
        zygosity = "GENO:GENO_0000458"
    if genotype == '1/1' or genotype == '1|1' or genotype == '1':
        zygosity = "GENO:GENO_0000136"

    return zygosity


""" Verifies if multiple samples are present, assign each genotype to sample"""


def multi_sample_case_level_data(case_level_data_biosampleId, data, row):
    # calculate number of samples
    samples = []
    genotypes = []
    zygozytes = []
    if "," in case_level_data_biosampleId:  # uses sample field to verify if multi sample
        samples = case_level_data_biosampleId.split(sep=",")  # split sample
    else:
        samples.append(case_level_data_biosampleId)
    if "," in case_level_data_biosampleId:  # uses sample field to verify if multi sample
        samples = case_level_data_biosampleId.split(sep=",")  # split sample
    else:
        samples.append(case_level_data_biosampleId)

    number_samples = len(samples)
    max_sample_position = 14 + number_samples
    case_level_data_zigosity_label_raw = row[14:max_sample_position]
    for genotype_raw in case_level_data_zigosity_label_raw:
        genotype = genotype_raw.split(sep=":")[0]
        genotypes.append(genotype)
        zygozyte = case_level_data(genotype)
        zygozytes.append(zygozyte)
    # add as many casa_level_data as we have samples
    for i in range(0, number_samples):
        sample = samples[i]
        genotype = genotypes[i]
        zygozyte = zygozytes[i]

        case_level_data_dict = {
            "zygosity": {
                "id": zygozyte,
                "label": genotype,
            },
            "biosampleId": sample
        }
        data["caseLevelData"].append(case_level_data_dict)
        # attempt to clear dictionary
        case_level_data_dict = {}

    return data


""" Receives Genotype, returns zygosity"""


def case_level_data(genotype):
    zygosity = None

    if genotype == '0/1' or genotype == '0|1' or genotype == '1/0' or genotype == '1|0' or genotype == '0':
        zygosity = "GENO:GENO_0000458"
    if genotype == '1/1' or genotype == '1|1' or genotype == '1':
        zygosity = "GENO:GENO_0000136"
    if genotype == "./." or genotype == ".|.":
        zygosity = None 
    # TODO "./."
    # unkown / unkown

    return zygosity


""" reads vcf """


def read_vcf():
    bcf_in = VariantFile("../data/delly.vcf.gz")  # auto-detect input format
    bcf_out = VariantFile('-', 'w', header=bcf_in.header)

    for rec in bcf_in.fetch('1', 10585, 14446418):
        bcf_out.write(rec)


""" Splits annotSV annotated TSV """


def pre_process_tsv():
    # counts number of samples
    samples_list = []
    with open('../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.tsv') as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        row1 = next(rd)
        samples = row1[7]
        samples_list = samples.split(sep=",")
        cut_end = len(samples_list)

    os.system(
        "cat ../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.tsv | cut -f1-29 | sort | uniq > ../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.uniq.tsv")


def get_value(row, position):
    value = row[position]
    return value


def read_tsv(data):
    with open('../data/lumpy.annotated.tsv') as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for row in rd:
            position_start_integer = get_value(row,2)
            position_end_integer = get_value(row,3)
            position_start = get_value(row,2)
            position_end = get_value(row,3)
            position_refseqId = get_value(row,1)
            position_assemblyId = "hg19"
            case_level_data_biosampleId = get_value(row,6)

            data = multi_sample_case_level_data(case_level_data_biosampleId, data, row)

            # calculate number of samples
            samples = []
            genotypes = []
            zygozytes = []
            if "," in case_level_data_biosampleId:  # uses sample field to verify if multi sample
                samples = case_level_data_biosampleId.split(sep=",")  # split sample
            else:
                samples.append(case_level_data_biosampleId)

            number_samples = len(samples)
            max_sample_position = 14 + number_samples
            case_level_data_zigosity_label_raw = row[14:max_sample_position]
            for genotype_raw in case_level_data_zigosity_label_raw:
                genotype = genotype_raw.split(sep=":")[0]
                genotypes.append(genotype)
                zygozyte = case_level_data(genotype)
                zygozytes.append(zygozyte)
            # add as many casa_level_data as we have samples
            for i in range(0, number_samples):
                sample = samples[i]
                genotype = genotypes[i]
                zygozyte = zygozytes[i]

                case_level_data_dict = {
                    "zygosity": {
                        "id": zygozyte,
                        "label": genotype,
                    },
                    "biosampleId": sample
                }
                print(case_level_data_dict)
                data["caseLevelData"].append(case_level_data_dict)
                # attempt to clear dictionary
                case_level_data_dict = {}

            variation_variant_type = get_value(row,5)
            variation_alternate_bases = get_value(row,9)
            variation_reference_bases = get_value(row,8)
            variation_location_interval_start_value = get_value(row,2)
            variation_location_interval_start_type = "Number"
            variation_location_interval_end_value = get_value(row,3)
            variation_location_interval_end_type = "Number"
            variation_location_type = "SequenceLocation"
            variation_location_sequence_id = get_value(row,5)
            variant_internal_id = get_value(row,0)
            identifiers_genomic_hgvs_id = get_value(row,5)
            variant_quality_filter = get_value(row,11)
            variant_quality_qual = get_value(row,10)
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
            filename = '../data/EEE_SV-Pop_1.ALL.sites.20181204.annotated.uniq.tsv'  # File name that was processed (vcf or tsv?)
            data['_info']['vcf2bff']['hostname'] = create_hostname()  # Hostname
            data['_info']['vcf2bff']['filein'] = create_filein(filename)  # File in
            data['_info']['vcf2bff']['user'] = create_user()  # User
            data['_info']['vcf2bff']['ncpuhost'] = create_ncpuhost()  # Number of cpu's
            data['_info']['vcf2bff']['fileout'] = create_fileout()  # File out
            data['_info']['vcf2bff']['cwd'] = create_cwd()  # Path to output folder
            data['_info']['vcf2bff']['cwd'] = create_cwd()  # Path to output folder
            data['_info']['vcf2bff']['projectDir'] = create_projectDir()  # Path to output final project directory
            data['_info']['vcf2bff']['version'] = create_version()  # Version of the Beacon
            write_json(data)
            data = read_json()


""" Adds semicomma to bff """


def bff_post_processing():
    os.system("sed -i 's/}{/},{/' ../results/g_variants_sv.json")


""" Annotates vcf with svAnnot """


def annotate_vcf(filename):
    ##add annotSV to path and run
    os.system(
        "export ANNOTSV=/home/mmoldes/Documents/EGA/bioteam/CNV_as_EGA_service/discovery_of_strucutral_variants_as_an_EGA_service/bin/AnnotSV && "
        "$ANNOTSV/bin/AnnotSV -SVinputFile /home/mmoldes/Documents/EGA/bioteam/structural_variants_beacon_v2/data/" + filename + " -genomeBuild GRCh37 -outputDir /home/mmoldes/Documents/EGA/bioteam/structural_variants_beacon_v2/data"
    )


""" write bff output in json array format """


def read_json():
    ## Transform from bff to python dictionary
    with open('../data/first_json_template.json') as json_file:
        data = json.load(json_file)  # Load json to python

    return data


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


def arguments_parser(args):
    parser = argparse.ArgumentParser(
        prog='sv_parser',
        description='What the program does',
        epilog='Text at the bottom of help')

    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max,
                        help='sum the integers (default: find the max)')

    return parser

    # TODO on main
    # args = parser.parse_args()
    print(args.filename, args.count, args.verbose)


if __name__ == "__main__":
    try:
        # configure logging
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s [in %(pathname)s:%(lineno)d]'
        logging.basicConfig(format=log_format)
        # execute main function
        # annotate_vcf("delly.vcf.gz") # ANNOTATES VCF WITH ANNOTSV
        # annotate_vcf("gridss.vcf.gz")
        # annotate_vcf("manta.vcf.gz")
        # annotate_vcf("lumpy.vcf.gz")
        # annotate_vcf("EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz")

        # todo
        # pre_process_tsv("manta.vcf.gz")

        data = read_json()  # READS BFF TEMPLATE
        read_tsv(data)  # POPULATE STUFF
        # internal_information(data)  # POPULATE INTERNAL INFORMATION
        bff_post_processing()

    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)
