#!/usr/bin/env python

"""sv_parse.py: Parses the annotsv annotation of a vcf file and generates a beacon compliant json array file ."""

__author__ = "Arnau Soler, Mauricio Moldes"
__version__ = "0.1"
__maintainer__ = "Arnau Soler, Mauricio Moldes"
__email__ = "arnau.soler@crg.eu, mauricio.moldes@crg.eu"
__status__ = "test"

import csv
import json
import math

from pysam import VariantFile
import logging
import sys
import os
import socket
import argparse
import pandas as pd

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


def parse_annotated_target_attribute_pandas(pd, attribute):
    try:
        value = pd[attribute]
        res = isNaN(value)
        if res:
            value = None
    except Exception as e:
        logger.error("Error: {}".format(e))
        value = None

    return value


def parse_annotated_target_attribute(annotated_tsv, collumn):
    try:
        value = annotated_tsv[collumn]
    except Exception as e:
        logger.error("Error: {}".format(e))
        value = None

    return value


""" PARSES ATTRIBUTE """


def parse_target_attribute(response_info, attribute):
    try:
        result = response_info[attribute]
    except Exception as e:
        logger.error("Error: {}".format(e))
        result = None
    return result


""" POPULATES DATA """


def parse_multiple_target_attribute(data, attribute, value):
    if attribute == "variant_internal_id":
        try:
            data['variantInternalId'] = value
        except Exception as e:
            logger.error("Error: {}".format(e))
    elif attribute == "identifiers_genomic_hgvs_id":
        try:
            data['identifiers']['genomicHGVSId'] = value
        except Exception as e:
            logger.error("Error: {}".format(e))
    elif attribute == "position_assemblyId":
        try:
            data['_position']['assemblyId'] = value
        except Exception as e:
            logger.error("Error: {}".format(e))
    elif attribute == "position_end":
        try:
            data['_position']['end'] = value
        except Exception as e:
            logger.error("Error: {}".format(e))

    return data


def parse_annotation():
    # specification documentation https://docs.genomebeacons.org/schemas-md/obj/molecularAttributes/

    return None


def read_pandas(filepath):
    # Read TSV file into DataFrame
    df = pd.read_table("../data/" + filepath + ".annotated.tsv")
    return df


# def ensg_glossary_traslator(term):
#     match term:
#         case 'Genome_annotation':
#             ensg_term = 'ENSGLOSSARY:0000001'
#         case 'Gene':
#             ensg_term = 'ENSGLOSSARY:0000002'
#         case 'Transcript':
#             ensg_term = 'ENSGLOSSARY:0000003'
#         case 'EST':
#             ensg_term = 'ENSGLOSSARY:0000004'
#         case 'Transcript_support_level':
#             ensg_term = 'ENSGLOSSARY:0000005'
#         case 'TSL_1':
#             ensg_term = 'ENSGLOSSARY:0000006'
#         case 'TSL_2':
#             ensg_term = 'ENSGLOSSARY:0000007'
#         case 'TSL_3':
#             ensg_term = 'ENSGLOSSARY:0000008'
#         case 'TSL_4':
#             ensg_term = 'ENSGLOSSARY:0000009'
#         case 'TSL_5':
#             ensg_term = 'ENSGLOSSARY:0000010'
#         case 'TSL_NA':
#             ensg_term = 'ENSGLOSSARY:0000011'
#         case 'APPRIS':
#             ensg_term = 'ENSGLOSSARY:0000012'
#         case 'APPRIS_P1':
#             ensg_term = 'ENSGLOSSARY:0000013'
#         case 'APPRIS_P2':
#             ensg_term = 'ENSGLOSSARY:0000014'
#         case 'APPRIS_P3':
#             ensg_term = 'ENSGLOSSARY:0000015'
#         case 'APPRIS_P4':
#             ensg_term = 'ENSGLOSSARY:0000016'
#         case 'APPRIS_P5':
#             ensg_term = 'ENSGLOSSARY:0000017'
#         case 'APPRIS_ALT1':
#             ensg_term = 'ENSGLOSSARY:0000018'
#         case 'APPRIS_ALT2':
#             ensg_term = 'ENSGLOSSARY:0000019'
#         case 'GENCODE_Basic':
#             ensg_term = 'ENSGLOSSARY:0000020'
#         case '5_incomplete' : ensg_term = 'ENSGLOSSARY:0000021'
#         case '3_incomplete' : ensg_term = 'ENSGLOSSARY:0000022'
#         case 'Ensembl_canonical':
#             ensg_term = 'ENSGLOSSARY:0000023'
#         case 'CCDS':
#             ensg_term = 'ENSGLOSSARY:0000024'
#         case 'Biotype':
#             ensg_term = 'ENSGLOSSARY:0000025'
#         case 'Protein_coding':
#             ensg_term = 'ENSGLOSSARY:0000026'
#         case 'Processed_transcript':
#             ensg_term = 'ENSGLOSSARY:0000027'
#         case 'Long_non-coding_RNA_(lncRNA)':
#             ensg_term = 'ENSGLOSSARY:0000028'
#         case 'Non_coding':
#             ensg_term = 'ENSGLOSSARY:0000029'
#         case '3_overlapping_ncRNA' : ensg_term = 'ENSGLOSSARY:0000030'
#         case 'Antisense':
#             ensg_term = 'ENSGLOSSARY:0000031'
#         case 'lincRNA_(long_intergenic_ncRNA)':
#             ensg_term = 'ENSGLOSSARY:0000032'
#         case 'Retained_intron':
#             ensg_term = 'ENSGLOSSARY:0000033'
#         case 'Sense_intronic':
#             ensg_term = 'ENSGLOSSARY:0000034'
#         case 'Sense_overlapping':
#             ensg_term = 'ENSGLOSSARY:0000035'
#         case 'Macro_lncRNA':
#             ensg_term = 'ENSGLOSSARY:0000036'
#         case 'ncRNA':
#             ensg_term = 'ENSGLOSSARY:0000037'
#         case 'miRNA':
#             ensg_term = 'ENSGLOSSARY:0000038'
#         case 'piRNA':
#             ensg_term = 'ENSGLOSSARY:0000039'
#         case 'rRNA':
#             ensg_term = 'ENSGLOSSARY:0000040'
#         case 'siRNA':
#             ensg_term = 'ENSGLOSSARY:0000041'
#         case 'snRNA':
#             ensg_term = 'ENSGLOSSARY:0000042'
#         case 'snoRNA':
#             ensg_term = 'ENSGLOSSARY:0000043'
#         case 'tRNA':
#             ensg_term = 'ENSGLOSSARY:0000044'
#         case 'vaultRNA':
#             ensg_term = 'ENSGLOSSARY:0000045'
#         case 'miscRNA':
#             ensg_term = 'ENSGLOSSARY:0000046'
#         case 'Pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000047'
#         case 'Processed_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000048'
#         case 'Unprocessed_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000049'
#         case 'Transcribed_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000050'
#         case 'Translated_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000051'
#         case 'Polymorphic_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000052'
#         case 'Unitary_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000053'
#         case 'IG_pseudogene':
#             ensg_term = 'ENSGLOSSARY:0000054'
#         case 'IG_gene':
#             ensg_term = 'ENSGLOSSARY:0000055'
#         case 'TR_gene':
#             ensg_term = 'ENSGLOSSARY:0000056'
#         case 'TEC_(To_be_Experimentally_Confirmed)':
#             ensg_term = 'ENSGLOSSARY:0000057'
#         case 'Readthrough':
#             ensg_term = 'ENSGLOSSARY:0000058'
#         case 'IG_V_gene':
#             ensg_term = 'ENSGLOSSARY:0000059'
#         case 'IG_D_gene':
#             ensg_term = 'ENSGLOSSARY:0000060'
#         case 'IG_J_gene':
#             ensg_term = 'ENSGLOSSARY:0000061'
#         case 'IG_C_gene':
#             ensg_term = 'ENSGLOSSARY:0000062'
#         case 'TR_V_gene':
#             ensg_term = 'ENSGLOSSARY:0000063'
#         case 'TR_D_gene':
#             ensg_term = 'ENSGLOSSARY:0000064'
#         case 'TR_J_gene':
#             ensg_term = 'ENSGLOSSARY:0000065'
#         case 'TR_C_gene':
#             ensg_term = 'ENSGLOSSARY:0000066'
#         case 'cDNA':
#             ensg_term = 'ENSGLOSSARY:0000067'
#         case 'CDS':
#             ensg_term = 'ENSGLOSSARY:0000068'
#         case 'Peptide':
#             ensg_term = 'ENSGLOSSARY:0000069'
#         case 'Protein_domain':
#             ensg_term = 'ENSGLOSSARY:0000070'
#         case 'Exon':
#             ensg_term = 'ENSGLOSSARY:0000071'
#         case 'Intron':
#             ensg_term = 'ENSGLOSSARY:0000072'
#         case 'Codon':
#             ensg_term = 'ENSGLOSSARY:0000073'
#         case 'Constitutive_exon':
#             ensg_term = 'ENSGLOSSARY:0000074'
#         case 'Phase':
#             ensg_term = 'ENSGLOSSARY:0000075'
#         case 'Flanking_sequence':
#             ensg_term = 'ENSGLOSSARY:0000076'
#         case 'Untranslated_region':
#             ensg_term = 'ENSGLOSSARY:0000077'
#         case '5_UTR' : ensg_term = 'ENSGLOSSARY:0000078'
#         case '3_UTR' : ensg_term = 'ENSGLOSSARY:0000079'
#         case 'Homologues':
#             ensg_term = 'ENSGLOSSARY:0000080'
#         case 'Gene_tree':
#             ensg_term = 'ENSGLOSSARY:0000081'
#         case 'Orthologues':
#             ensg_term = 'ENSGLOSSARY:0000082'
#         case '1-to-1_orthologues':
#             ensg_term = 'ENSGLOSSARY:0000083'
#         case '1-to-many_orthologues':
#             ensg_term = 'ENSGLOSSARY:0000084'
#         case 'Many-to-many_orthologues':
#             ensg_term = 'ENSGLOSSARY:0000085'
#         case 'Paralogues':
#             ensg_term = 'ENSGLOSSARY:0000086'
#         case 'Between_species_paralogues':
#             ensg_term = 'ENSGLOSSARY:0000087'
#         case 'Gene_split':
#             ensg_term = 'ENSGLOSSARY:0000088'
#         case 'Other_paralogues':
#             ensg_term = 'ENSGLOSSARY:0000089'
#         case 'Within_species_paralogues':
#             ensg_term = 'ENSGLOSSARY:0000090'
#         case 'Homoeologues':
#             ensg_term = 'ENSGLOSSARY:0000091'
#         case 'Variant':
#             ensg_term = 'ENSGLOSSARY:0000092'
#         case 'QTL':
#             ensg_term = 'ENSGLOSSARY:0000093'
#         case 'eQTL':
#             ensg_term = 'ENSGLOSSARY:0000094'
#         case 'Evidence_status':
#             ensg_term = 'ENSGLOSSARY:0000095'
#         case 'Sequence_variant':
#             ensg_term = 'ENSGLOSSARY:0000096'
#         case 'Structural_variant':
#             ensg_term = 'ENSGLOSSARY:0000097'
#         case 'SNP':
#             ensg_term = 'ENSGLOSSARY:0000098'
#         case 'Insertion':
#             ensg_term = 'ENSGLOSSARY:0000099'
#         case 'Deletion':
#             ensg_term = 'ENSGLOSSARY:0000100'
#         case 'Indel':
#             ensg_term = 'ENSGLOSSARY:0000101'
#         case 'Substitution':
#             ensg_term = 'ENSGLOSSARY:0000102'
#         case 'CNV':
#             ensg_term = 'ENSGLOSSARY:0000103'
#         case 'Inversion':
#             ensg_term = 'ENSGLOSSARY:0000104'
#         case 'Translocation':
#             ensg_term = 'ENSGLOSSARY:0000105'
#         case 'Allele_(variant)':
#             ensg_term = 'ENSGLOSSARY:0000106'
#         case 'Allele_(gene)':
#             ensg_term = 'ENSGLOSSARY:0000107'
#         case 'Reference_allele':
#             ensg_term = 'ENSGLOSSARY:0000108'
#         case 'Alternative_allele':
#             ensg_term = 'ENSGLOSSARY:0000109'
#         case 'Major_allele':
#             ensg_term = 'ENSGLOSSARY:0000110'
#         case 'Minor_allele':
#             ensg_term = 'ENSGLOSSARY:0000111'
#         case 'Private_allele':
#             ensg_term = 'ENSGLOSSARY:0000112'
#         case 'Ancestral_allele':
#             ensg_term = 'ENSGLOSSARY:0000113'
#         case 'Minor_allele_frequency':
#             ensg_term = 'ENSGLOSSARY:0000114'
#         case 'Highest_population_MAF':
#             ensg_term = 'ENSGLOSSARY:0000115'
#         case 'Global_MAF':
#             ensg_term = 'ENSGLOSSARY:0000116'
#         case 'Genotype':
#             ensg_term = 'ENSGLOSSARY:0000117'
#         case 'Genetic_marker':
#             ensg_term = 'ENSGLOSSARY:0000118'
#         case 'Tandem_repeat':
#             ensg_term = 'ENSGLOSSARY:0000119'
#         case 'Alu_insertion':
#             ensg_term = 'ENSGLOSSARY:0000120'
#         case 'Complex_structural_alteration':
#             ensg_term = 'ENSGLOSSARY:0000121'
#         case 'Complex_substitution':
#             ensg_term = 'ENSGLOSSARY:0000122'
#         case 'Interchromosomal_breakpoint':
#             ensg_term = 'ENSGLOSSARY:0000123'
#         case 'Interchromosomal_translocation':
#             ensg_term = 'ENSGLOSSARY:0000124'
#         case 'Intrachromosomal_breakpoint':
#             ensg_term = 'ENSGLOSSARY:0000125'
#         case 'Intrachromosomal_translocation':
#             ensg_term = 'ENSGLOSSARY:0000126'
#         case 'Loss_of_heterozygosity':
#             ensg_term = 'ENSGLOSSARY:0000127'
#         case 'Mobile_element_deletion':
#             ensg_term = 'ENSGLOSSARY:0000128'
#         case 'Mobile_element_insertion':
#             ensg_term = 'ENSGLOSSARY:0000129'
#         case 'Novel_sequence_insertion':
#             ensg_term = 'ENSGLOSSARY:0000130'
#         case 'Short_tandem_repeat_variant':
#             ensg_term = 'ENSGLOSSARY:0000131'
#         case 'Tandem_duplication':
#             ensg_term = 'ENSGLOSSARY:0000132'
#         case 'Probe':
#             ensg_term = 'ENSGLOSSARY:0000133'
#         case 'Variant_consequence':
#             ensg_term = 'ENSGLOSSARY:0000134'
#         case 'Variant_impact':
#             ensg_term = 'ENSGLOSSARY:0000135'
#         case 'High_impact_variant_consequence':
#             ensg_term = 'ENSGLOSSARY:0000136'
#         case 'Moderate_impact_variant_consequence':
#             ensg_term = 'ENSGLOSSARY:0000137'
#         case 'Low_impact_variant_consequence':
#             ensg_term = 'ENSGLOSSARY:0000138'
#         case 'Modifier_impact_variant_consequence':
#             ensg_term = 'ENSGLOSSARY:0000139'
#         case 'Transcript_ablation':
#             ensg_term = 'ENSGLOSSARY:0000140'
#         case 'Splice_acceptor_variant':
#             ensg_term = 'ENSGLOSSARY:0000141'
#         case 'Splice_donor_variant':
#             ensg_term = 'ENSGLOSSARY:0000142'
#         case 'Stop_gained':
#             ensg_term = 'ENSGLOSSARY:0000143'
#         case 'Frameshift_variant':
#             ensg_term = 'ENSGLOSSARY:0000144'
#         case 'Stop_lost':
#             ensg_term = 'ENSGLOSSARY:0000145'
#         case 'Start_lost':
#             ensg_term = 'ENSGLOSSARY:0000146'
#         case 'Transcript_amplification':
#             ensg_term = 'ENSGLOSSARY:0000147'
#         case 'Inframe_insertion':
#             ensg_term = 'ENSGLOSSARY:0000148'
#         case 'Inframe_deletion':
#             ensg_term = 'ENSGLOSSARY:0000149'
#         case 'Missense_variant':
#             ensg_term = 'ENSGLOSSARY:0000150'
#         case 'Protein_altering_variant':
#             ensg_term = 'ENSGLOSSARY:0000151'
#         case 'Splice_region_variant':
#             ensg_term = 'ENSGLOSSARY:0000152'
#         case 'Incomplete_terminal_codon_variant':
#             ensg_term = 'ENSGLOSSARY:0000153'
#         case 'Stop_retained_variant':
#             ensg_term = 'ENSGLOSSARY:0000154'
#         case 'Synonymous_variant':
#             ensg_term = 'ENSGLOSSARY:0000155'
#         case 'Coding_sequence_variant':
#             ensg_term = 'ENSGLOSSARY:0000156'
#         case 'Mature_miRNA_variant':
#             ensg_term = 'ENSGLOSSARY:0000157'
#         case '5_prime_UTR_variant':
#             ensg_term = 'ENSGLOSSARY:0000158'
#         case '3_prime_UTR_variant':
#             ensg_term = 'ENSGLOSSARY:0000159'
#         case 'Non_coding_transcript_exon_variant':
#             ensg_term = 'ENSGLOSSARY:0000160'
#         case 'Intron_variant':
#             ensg_term = 'ENSGLOSSARY:0000161'
#         case 'NMD_transcript_variant':
#             ensg_term = 'ENSGLOSSARY:0000162'
#         case 'Non_coding_transcript_variant':
#             ensg_term = 'ENSGLOSSARY:0000163'
#         case 'Upstream_gene_variant':
#             ensg_term = 'ENSGLOSSARY:0000164'
#         case 'Downstream_gene_variant':
#             ensg_term = 'ENSGLOSSARY:0000165'
#         case 'TFBS_ablation':
#             ensg_term = 'ENSGLOSSARY:0000166'
#         case 'TFBS_amplification':
#             ensg_term = 'ENSGLOSSARY:0000167'
#         case 'TF_binding_site_variant':
#             ensg_term = 'ENSGLOSSARY:0000168'
#         case 'Regulatory_region_ablation':
#             ensg_term = 'ENSGLOSSARY:0000169'
#         case 'Regulatory_region_amplification':
#             ensg_term = 'ENSGLOSSARY:0000170'
#         case 'Feature_elongation':
#             ensg_term = 'ENSGLOSSARY:0000171'
#         case 'Regulatory_region_variant':
#             ensg_term = 'ENSGLOSSARY:0000172'
#         case 'Feature_truncation':
#             ensg_term = 'ENSGLOSSARY:0000173'
#         case 'Intergenic_variant':
#             ensg_term = 'ENSGLOSSARY:0000174'
#         case 'Ambiguity_code':
#             ensg_term = 'ENSGLOSSARY:0000175'
#         case 'Flagged_variant':
#             ensg_term = 'ENSGLOSSARY:0000176'
#         case 'Clinical_significance':
#             ensg_term = 'ENSGLOSSARY:0000177'
#         case 'Linkage_disequilibrium':
#             ensg_term = 'ENSGLOSSARY:0000178'
#         case 'r2':
#             ensg_term = 'ENSGLOSSARY:0000179'
#         case 'D' : ensg_term = 'ENSGLOSSARY:0000180'
#         case 'Haplotype_(variation)':
#             ensg_term = 'ENSGLOSSARY:0000181'
#         case 'Transcript_haplotype':
#             ensg_term = 'ENSGLOSSARY:0000182'
#         case 'Genome':
#             ensg_term = 'ENSGLOSSARY:0000183'
#         case 'Genome_assembly':
#             ensg_term = 'ENSGLOSSARY:0000184'
#         case 'Coverage':
#             ensg_term = 'ENSGLOSSARY:0000185'
#         case 'Primary_assembly':
#             ensg_term = 'ENSGLOSSARY:0000186'
#         case 'Alternative_sequence':
#             ensg_term = 'ENSGLOSSARY:0000187'
#         case 'Patch':
#             ensg_term = 'ENSGLOSSARY:0000188'
#         case 'Haplotype_(genome)':
#             ensg_term = 'ENSGLOSSARY:0000189'
#         case 'Novel_patch':
#             ensg_term = 'ENSGLOSSARY:0000190'
#         case 'Fix_patch':
#             ensg_term = 'ENSGLOSSARY:0000191'
#         case 'Contig':
#             ensg_term = 'ENSGLOSSARY:0000192'
#         case 'Scaffold':
#             ensg_term = 'ENSGLOSSARY:0000193'
#         case 'Cytogenetic_band':
#             ensg_term = 'ENSGLOSSARY:0000194'
#         case 'Clone':
#             ensg_term = 'ENSGLOSSARY:0000195'
#         case 'BAC':
#             ensg_term = 'ENSGLOSSARY:0000196'
#         case 'YAC':
#             ensg_term = 'ENSGLOSSARY:0000197'
#         case 'Cosmid':
#             ensg_term = 'ENSGLOSSARY:0000198'
#         case 'Base_pairs_(genome_size)':
#             ensg_term = 'ENSGLOSSARY:0000199'
#         case 'Golden_path_(genome_size)':
#             ensg_term = 'ENSGLOSSARY:0000200'
#         case 'Coordinate_system':
#             ensg_term = 'ENSGLOSSARY:0000201'
#         case 'Karyotype':
#             ensg_term = 'ENSGLOSSARY:0000202'
#         case 'PAR':
#             ensg_term = 'ENSGLOSSARY:0000203'
#         case 'Slice':
#             ensg_term = 'ENSGLOSSARY:0000204'
#         case 'Toplevel':
#             ensg_term = 'ENSGLOSSARY:0000205'
#         case 'Placed_scaffold':
#             ensg_term = 'ENSGLOSSARY:0000206'
#         case 'Unplaced_scaffold':
#             ensg_term = 'ENSGLOSSARY:0000207'
#         case 'Ensembl_sources':
#             ensg_term = 'ENSGLOSSARY:0000208'
#         case 'Gene_source_database':
#             ensg_term = 'ENSGLOSSARY:0000209'
#         case 'GENCODE':
#             ensg_term = 'ENSGLOSSARY:0000210'
#         case 'RefSeq':
#             ensg_term = 'ENSGLOSSARY:0000211'
#         case 'UniProt':
#             ensg_term = 'ENSGLOSSARY:0000212'
#         case 'IMGT':
#             ensg_term = 'ENSGLOSSARY:0000213'
#         case 'SwissProt':
#             ensg_term = 'ENSGLOSSARY:0000214'
#         case 'TrEMBL':
#             ensg_term = 'ENSGLOSSARY:0000215'
#         case 'INSDC':
#             ensg_term = 'ENSGLOSSARY:0000216'
#         case 'ENA':
#             ensg_term = 'ENSGLOSSARY:0000217'
#         case 'GenBank_(database)':
#             ensg_term = 'ENSGLOSSARY:0000218'
#         case 'DDBJ':
#             ensg_term = 'ENSGLOSSARY:0000219'
#         case 'Gene_Ontology':
#             ensg_term = 'ENSGLOSSARY:0000220'
#         case 'HGNC':
#             ensg_term = 'ENSGLOSSARY:0000221'
#         case 'Rfam':
#             ensg_term = 'ENSGLOSSARY:0000222'
#         case 'miRbase':
#             ensg_term = 'ENSGLOSSARY:0000223'
#         case 'MGI':
#             ensg_term = 'ENSGLOSSARY:0000224'
#         case 'zFIN':
#             ensg_term = 'ENSGLOSSARY:0000225'
#         case 'SGD':
#             ensg_term = 'ENSGLOSSARY:0000226'
#         case 'UCSC_Genome_Browser':
#             ensg_term = 'ENSGLOSSARY:0000227'
#         case 'Epigenome_source_database':
#             ensg_term = 'ENSGLOSSARY:0000228'
#         case 'ENCODE':
#             ensg_term = 'ENSGLOSSARY:0000229'
#         case 'Blueprint_Epigenomes':
#             ensg_term = 'ENSGLOSSARY:0000230'
#         case 'Roadmap_Epigenomics':
#             ensg_term = 'ENSGLOSSARY:0000231'
#         case 'Variation_source_database':
#             ensg_term = 'ENSGLOSSARY:0000232'
#         case 'dbSNP':
#             ensg_term = 'ENSGLOSSARY:0000233'
#         case 'EVA':
#             ensg_term = 'ENSGLOSSARY:0000234'
#         case 'dbVar':
#             ensg_term = 'ENSGLOSSARY:0000235'
#         case 'DGVa':
#             ensg_term = 'ENSGLOSSARY:0000236'
#         case '1000_Genomes_project':
#             ensg_term = 'ENSGLOSSARY:0000237'
#         case 'gnomAD':
#             ensg_term = 'ENSGLOSSARY:0000238'
#         case 'TOPMed':
#             ensg_term = 'ENSGLOSSARY:0000239'
#         case 'UK10K':
#             ensg_term = 'ENSGLOSSARY:0000240'
#         case 'HapMap':
#             ensg_term = 'ENSGLOSSARY:0000241'
#         case 'ClinVar':
#             ensg_term = 'ENSGLOSSARY:0000242'
#         case 'Phenotype_source_database':
#             ensg_term = 'ENSGLOSSARY:0000243'
#         case 'OMIM':
#             ensg_term = 'ENSGLOSSARY:0000244'
#         case 'OMIA':
#             ensg_term = 'ENSGLOSSARY:0000245'
#         case 'Orphanet':
#             ensg_term = 'ENSGLOSSARY:0000246'
#         case 'GWAS_catalog':
#             ensg_term = 'ENSGLOSSARY:0000247'
#         case 'IMPC':
#             ensg_term = 'ENSGLOSSARY:0000248'
#         case 'Animal_QTLdb':
#             ensg_term = 'ENSGLOSSARY:0000249'
#         case 'HGMD':
#             ensg_term = 'ENSGLOSSARY:0000250'
#         case 'COSMIC':
#             ensg_term = 'ENSGLOSSARY:0000251'
#         case 'Protein_source_database':
#             ensg_term = 'ENSGLOSSARY:0000252'
#         case 'PDB':
#             ensg_term = 'ENSGLOSSARY:0000253'
#         case 'Algorithm':
#             ensg_term = 'ENSGLOSSARY:0000254'
#         case 'Ensembl_Genebuild':
#             ensg_term = 'ENSGLOSSARY:0000255'
#         case 'Ensembl_Havana':
#             ensg_term = 'ENSGLOSSARY:0000256'
#         case 'Ensembl_Regulatory_Build':
#             ensg_term = 'ENSGLOSSARY:0000257'
#         case 'Ensembl_gene_tree_pipeline':
#             ensg_term = 'ENSGLOSSARY:0000258'
#         case 'InterProScan':
#             ensg_term = 'ENSGLOSSARY:0000259'
#         case 'BLAST':
#             ensg_term = 'ENSGLOSSARY:0000260'
#         case 'BLAT':
#             ensg_term = 'ENSGLOSSARY:0000261'
#         case 'DUST':
#             ensg_term = 'ENSGLOSSARY:0000262'
#         case 'Eponine':
#             ensg_term = 'ENSGLOSSARY:0000263'
#         case 'GeneWise':
#             ensg_term = 'ENSGLOSSARY:0000264'
#         case 'Exonerate':
#             ensg_term = 'ENSGLOSSARY:0000265'
#         case 'Projection_build':
#             ensg_term = 'ENSGLOSSARY:0000266'
#         case 'GENSCAN':
#             ensg_term = 'ENSGLOSSARY:0000267'
#         case 'SIFT':
#             ensg_term = 'ENSGLOSSARY:0000268'
#         case 'PolyPhen':
#             ensg_term = 'ENSGLOSSARY:0000269'
#         case 'RepeatMasker':
#             ensg_term = 'ENSGLOSSARY:0000270'
#         case 'BLOSUM_62':
#             ensg_term = 'ENSGLOSSARY:0000271'
#         case 'VEP':
#             ensg_term = 'ENSGLOSSARY:0000272'
#         case 'File_formats':
#             ensg_term = 'ENSGLOSSARY:0000273'
#         case 'HGVS_nomenclature':
#             ensg_term = 'ENSGLOSSARY:0000274'
#         case 'VCF':
#             ensg_term = 'ENSGLOSSARY:0000275'
#         case 'BED':
#             ensg_term = 'ENSGLOSSARY:0000276'
#         case 'FASTA':
#             ensg_term = 'ENSGLOSSARY:0000277'
#         case 'BAM/CRAM':
#             ensg_term = 'ENSGLOSSARY:0000278'
#         case 'BigBed':
#             ensg_term = 'ENSGLOSSARY:0000279'
#         case 'Ensembl_default_(VEP)':
#             ensg_term = 'ENSGLOSSARY:0000280'
#         case 'BedGraph':
#             ensg_term = 'ENSGLOSSARY:0000281'
#         case 'GTF':
#             ensg_term = 'ENSGLOSSARY:0000282'
#         case 'GFF':
#             ensg_term = 'ENSGLOSSARY:0000283'
#         case 'PSL':
#             ensg_term = 'ENSGLOSSARY:0000284'
#         case 'Wiggle':
#             ensg_term = 'ENSGLOSSARY:0000285'
#         case 'BigWig':
#             ensg_term = 'ENSGLOSSARY:0000286'
#         case 'Pairwise_interactions_(WashU)':
#             ensg_term = 'ENSGLOSSARY:0000287'
#         case 'chain':
#             ensg_term = 'ENSGLOSSARY:0000288'
#         case 'Newick':
#             ensg_term = 'ENSGLOSSARY:0000289'
#         case 'EMBL_(file_format)':
#             ensg_term = 'ENSGLOSSARY:0000290'
#         case 'GenBank_(file_format)':
#             ensg_term = 'ENSGLOSSARY:0000291'
#         case 'EMF_Alignment_format':
#             ensg_term = 'ENSGLOSSARY:0000292'
#         case 'MAF':
#             ensg_term = 'ENSGLOSSARY:0000293'
#         case 'MySQL':
#             ensg_term = 'ENSGLOSSARY:0000294'
#         case 'VEP_cache':
#             ensg_term = 'ENSGLOSSARY:0000295'
#         case 'GVF':
#             ensg_term = 'ENSGLOSSARY:0000296'
#         case 'PhyloXML':
#             ensg_term = 'ENSGLOSSARY:0000297'
#         case 'OrthoXML':
#             ensg_term = 'ENSGLOSSARY:0000298'
#         case 'RDF':
#             ensg_term = 'ENSGLOSSARY:0000299'
#         case 'AGP':
#             ensg_term = 'ENSGLOSSARY:0000300'
#         case 'Repeat':
#             ensg_term = 'ENSGLOSSARY:0000301'
#         case 'Repeat_masking':
#             ensg_term = 'ENSGLOSSARY:0000302'
#         case 'Hard_masked':
#             ensg_term = 'ENSGLOSSARY:0000303'
#         case 'Soft_masked':
#             ensg_term = 'ENSGLOSSARY:0000304'
#         case 'Alu_insertion':
#             ensg_term = 'ENSGLOSSARY:0000305'
#         case 'Microsatellite':
#             ensg_term = 'ENSGLOSSARY:0000306'
#         case 'Centromere':
#             ensg_term = 'ENSGLOSSARY:0000307'
#         case 'Low_complexity_regions':
#             ensg_term = 'ENSGLOSSARY:0000308'
#         case 'RNA_repeats':
#             ensg_term = 'ENSGLOSSARY:0000309'
#         case 'Satellite_repeats':
#             ensg_term = 'ENSGLOSSARY:0000310'
#         case 'Simple_repeats':
#             ensg_term = 'ENSGLOSSARY:0000311'
#         case 'Tandem_repeats':
#             ensg_term = 'ENSGLOSSARY:0000312'
#         case 'LTRs':
#             ensg_term = 'ENSGLOSSARY:0000313'
#         case 'Type_I_Transposons/LINE':
#             ensg_term = 'ENSGLOSSARY:0000314'
#         case 'Type_I_Transposons/SINE':
#             ensg_term = 'ENSGLOSSARY:0000315'
#         case 'Type_II_Transposons':
#             ensg_term = 'ENSGLOSSARY:0000316'
#         case 'Unknown_repeat':
#             ensg_term = 'ENSGLOSSARY:0000317'
#         case 'Alignments':
#             ensg_term = 'ENSGLOSSARY:0000318'
#         case 'Whole_genome_alignment':
#             ensg_term = 'ENSGLOSSARY:0000319'
#         case 'Pairwise_whole_genome_alignment':
#             ensg_term = 'ENSGLOSSARY:0000320'
#         case 'Multiple_whole_genome_alignment':
#             ensg_term = 'ENSGLOSSARY:0000321'
#         case 'Synteny':
#             ensg_term = 'ENSGLOSSARY:0000322'
#         case 'CIGAR':
#             ensg_term = 'ENSGLOSSARY:0000323'
#         case 'Identity':
#             ensg_term = 'ENSGLOSSARY:0000324'
#         case 'Wasabi':
#             ensg_term = 'ENSGLOSSARY:0000325'
#         case 'Similarity':
#             ensg_term = 'ENSGLOSSARY:0000326'
#         case 'Pecan':
#             ensg_term = 'ENSGLOSSARY:0000327'
#         case 'EPO':
#             ensg_term = 'ENSGLOSSARY:0000328'
#         case 'Progressive_cactus':
#             ensg_term = 'ENSGLOSSARY:0000329'
#         case 'LastZ':
#             ensg_term = 'ENSGLOSSARY:0000330'
#         case 'BlastZ':
#             ensg_term = 'ENSGLOSSARY:0000331'
#         case 'Translated_Blat':
#             ensg_term = 'ENSGLOSSARY:0000332'
#         case 'Regulatory_features':
#             ensg_term = 'ENSGLOSSARY:0000333'
#         case 'Promoters':
#             ensg_term = 'ENSGLOSSARY:0000334'
#         case 'Promoter_flanking_regions':
#             ensg_term = 'ENSGLOSSARY:0000335'
#         case 'Enhancers':
#             ensg_term = 'ENSGLOSSARY:0000336'
#         case 'CTCF_binding_sites':
#             ensg_term = 'ENSGLOSSARY:0000337'
#         case 'Transcription_factor_binding_sites':
#             ensg_term = 'ENSGLOSSARY:0000338'
#         case 'Open_chromatin_regions':
#             ensg_term = 'ENSGLOSSARY:0000339'
#         case 'Regulatory_activity':
#             ensg_term = 'ENSGLOSSARY:0000340'
#         case 'Active':
#             ensg_term = 'ENSGLOSSARY:0000341'
#         case 'Poised':
#             ensg_term = 'ENSGLOSSARY:0000342'
#         case 'Repressed':
#             ensg_term = 'ENSGLOSSARY:0000343'
#         case 'Inactive':
#             ensg_term = 'ENSGLOSSARY:0000344'
#         case 'NA':
#             ensg_term = 'ENSGLOSSARY:0000345'
#         case 'Epigenome_evidence':
#             ensg_term = 'ENSGLOSSARY:0000346'
#         case 'ChIP-seq':
#             ensg_term = 'ENSGLOSSARY:0000347'
#         case 'DNase_sensitivity':
#             ensg_term = 'ENSGLOSSARY:0000348'
#         case 'Transcription_factor':
#             ensg_term = 'ENSGLOSSARY:0000349'
#         case 'Histone_modification':
#             ensg_term = 'ENSGLOSSARY:0000350'
#         case 'DNA_methylation':
#             ensg_term = 'ENSGLOSSARY:0000351'
#         case 'Bisulfite_sequencing':
#             ensg_term = 'ENSGLOSSARY:0000352'
#         case 'Signal':
#             ensg_term = 'ENSGLOSSARY:0000353'
#         case 'Peak':
#             ensg_term = 'ENSGLOSSARY:0000354'
#         case 'Transcription_factor_binding_motif':
#             ensg_term = 'ENSGLOSSARY:0000355'
#         case 'Epigenome':
#             ensg_term = 'ENSGLOSSARY:0000356'
#         case 'Marker':
#             ensg_term = 'ENSGLOSSARY:0000357'
#         case 'UniSTS':
#             ensg_term = 'ENSGLOSSARY:0000358'
#         case 'External_reference':
#             ensg_term = 'ENSGLOSSARY:0000359'
#         case 'CADD':
#             ensg_term = 'ENSGLOSSARY:0000360'
#         case 'REVEL':
#             ensg_term = 'ENSGLOSSARY:0000361'
#         case 'MutationAssessor':
#             ensg_term = 'ENSGLOSSARY:0000362'
#         case 'MetaLR':
#             ensg_term = 'ENSGLOSSARY:0000363'
#         case 'MANE':
#             ensg_term = 'ENSGLOSSARY:0000364'
#         case 'MANE_Select':
#             ensg_term = 'ENSGLOSSARY:0000365'
#         case 'TAGENE':
#             ensg_term = 'ENSGLOSSARY:0000367'
#         case 'Stop_codon_readthrough':
#             ensg_term = 'ENSGLOSSARY:0000368'
#         case 'Forward_strand':
#             ensg_term = 'ENSGLOSSARY:0000369'
#         case 'Reverse_strand':
#             ensg_term = 'ENSGLOSSARY:0000370'
#         case 'RefSeq_Match':
#             ensg_term = 'ENSGLOSSARY:0000371'
#         case 'UniProt_Match':
#             ensg_term = 'ENSGLOSSARY:0000372'
#         case 'Nonsense_Mediated_Decay':
#             ensg_term = 'ENSGLOSSARY:0000373'
#         case 'Non-ATG_start':
#             ensg_term = 'ENSGLOSSARY:0000374'
#         case 'MANE_Plus_Clinical':
#             ensg_term = 'ENSGLOSSARY:0000375'
#         case 'GENCODE_Comprehensive':
#             ensg_term = 'ENSGLOSSARY:0000376'
#         case _:
#             ensg_term = "Not parsable"
#     return ensg_term


def isNaN(string):
    return string != string


def parse_pandas(pd, data, filepath):
    for index, row in pd.iterrows():
        print(row['AnnotSV_ID'], row['SV_chrom'])  # gives a sense of progression
        ronw_annot_sv_id = row['AnnotSV_ID']
        if ronw_annot_sv_id == "1_91757755_91767193_DUP_1":
            print("breakpoint!")
        # debug !

        position_start_integer = parse_annotated_target_attribute_pandas(row, 'SV_start')
        position_end_integer = parse_annotated_target_attribute_pandas(row, 'SV_end')
        position_start = parse_annotated_target_attribute_pandas(row, 'SV_start')
        position_end = parse_annotated_target_attribute_pandas(row, 'SV_end')
        position_refseqId = parse_annotated_target_attribute_pandas(row, 'SV_chrom')
        position_assemblyId = "hg19"
        case_level_data_biosampleId = parse_annotated_target_attribute_pandas(row, 'Samples_ID')

        genes = parse_annotated_target_attribute_pandas(row, 'Gene_name')

        if genes is None:
            # if math.isnan(genes):
            data['molecularAttributes']['geneIds'][0] = None
        else:
            list_genes = []
            if ";" in genes:
                list_genes = genes.split(";")
            else:
                list_genes.append(genes)

            set_of_genes = ','.join(list_genes)

            data['molecularAttributes']['geneIds'][0] = set_of_genes

        # calculate number of samples
        samples = []
        genotypes = []
        zygozytes = []
        if "," in case_level_data_biosampleId:  # uses sample field to verify if multi sample
            samples = case_level_data_biosampleId.split(sep=",")  # split sample
        else:
            samples.append(case_level_data_biosampleId)

        for sample in samples:
            case_level_data_zigosity_label_raw = row[sample]

            genotype = case_level_data_zigosity_label_raw.split(sep=":")[0]
            genotypes.append(genotype)
            zygozyte = case_level_data(genotype)
            zygozytes.append(zygozyte)

        number_samples = len(samples)

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

        variation_variant_type = parse_annotated_target_attribute_pandas(row, 'SV_type')
        variation_alternate_bases = parse_annotated_target_attribute_pandas(row, 'ALT')
        variation_reference_bases = parse_annotated_target_attribute_pandas(row, 'REF')
        variation_location_interval_start_value = parse_annotated_target_attribute_pandas(row, 'SV_start')
        variation_location_interval_start_type = "Number"
        variation_location_interval_end_value = parse_annotated_target_attribute_pandas(row, 'SV_end')
        variation_location_interval_end_type = "Number"
        variation_location_type = "SequenceLocation"
        variation_location_sequence_id = parse_annotated_target_attribute_pandas(row, 'ID')
        variant_internal_id = parse_annotated_target_attribute_pandas(row, 'AnnotSV_ID')
        identifiers_genomic_hgvs_id = parse_annotated_target_attribute_pandas(row, 5)
        variant_quality_filter = parse_annotated_target_attribute_pandas(row, 'FILTER')
        variant_quality_qual = parse_annotated_target_attribute_pandas(row, 'QUAL')
        variation_location_interval_type = "SequenceLocation"

        AnnotSV_ranking_score = parse_annotated_target_attribute_pandas(row, 'AnnotSV_ranking_score')
        AnnotSV_ranking_criteria = parse_annotated_target_attribute_pandas(row, 'AnnotSV_ranking_criteria')
        ACMG_class = parse_annotated_target_attribute_pandas(row, 'ACMG_class')
        converted_acmg_class = convert_acmg_class(ACMG_class)

        data['molecularAttributes']['annotationImpact'] = converted_acmg_class
        # data['identifiers']['genomicHGVSId'] = identifiers_genomic_hgvs_id
        # data['_position']['assemblyId'] = position_assemblyId

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
        write_json(data, filepath)
        data = read_json()


""" Adds semicomma to bff """


def bff_post_processing(filename):
    os.system("sed -i 's/}{/},{/' ../results/" + filename + ".g_variants_sv.json")


""" Annotates vcf with svAnnot """


def annotate_vcf(filename, genome_build, annotation_mode):
    ##add annotSV to path and run
    os.system(
        "export ANNOTSV=/home/mmoldes/Documents/EGA/bioteam/CNV_as_EGA_service/discovery_of_strucutral_variants_as_an_EGA_service/bin/AnnotSV && "
        "$ANNOTSV/bin/AnnotSV -SVinputFile /home/mmoldes/Documents/EGA/bioteam/structural_variants_beacon_v2/data/" + filename + ".vcf.gz -genomeBuild " + genome_build + " -outputDir /home/mmoldes/Documents/EGA/bioteam/structural_variants_beacon_v2/data -annotationMode " + annotation_mode + "")


""" write bff output in json array format """


def read_json():
    ## Transform from bff to python dictionary
    with open('../data/first_json_template.json') as json_file:
        data = json.load(json_file)  # Load json to python

    return data


""" write bff output in json array format """


def write_json(data, filepath):
    # Serializing json
    json_object = json.dumps(data, indent=4)

    # Writing to sample.json
    with open("../results/" + filepath + ".g_variants_sv.json", "a") as outfile:
        outfile.write(json_object)


def convert_acmg_class(acmg_value):
    """
    SV ranking class into 1 of 5:
    class 1 (benign)
    class 2 (likely benign)
    class 3 (variant of unknown significance)
    class 4 (likely pathogenic)
    class 5 (pathogenic)
    class NA (Non Attributed)
    """
    if acmg_value is None:
        result = None
    elif acmg_value == 1.0:
        result = "benign"
    elif acmg_value == 2.0:
        result = "likely benign"
    elif acmg_value == 3.0:
        result = "variant of unknown significance"
    elif acmg_value == 4.0:
        result = "likely pathogenic"
    elif acmg_value == 5.0:
        result = "pathogenic"
    elif acmg_value == "NA":
        result = "Non Attributed"

    return result


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

        # filenames = ["nstd137.GRCh37.variant_call",
        #             "nstd152.GRCh37.variant_call", "nstd162.GRCh37.variant_call",
        #             "nstd167.GRCh37.variant_call", "nstd171.GRCh37.variant_call", "nstd175.GRCh37.variant_call",
        #             "nstd186.GRCh37.variant_call", "nstd102.GRCh37.variant_call"]

        # filenames = ["HG001_GRCh37_1_22_v4.2.1_benchmark",
        #              "HG002_GRCh37_1_22_v4.2.1_benchmark",
        #              "HG003_GRCh37_1_22_v4.2.1_benchmark",
        #              "HG004_GRCh37_1_22_v4.2.1_benchmark",
        #              "HG005_GRCh37_1_22_v4.2.1_benchmark",
        #              "HG007_GRCh37_1_22_v4.2.1_benchmark"]

        filenames = ["HG007_GRCh37_1_22_v4.2.1_benchmark"]

        # filenames = [            "all_ins_raw",            "all_trp_raw",            "all_tra_raw",            "all_middel_raw",            "all_inv_raw"        ]

        # filenames = ["gridss","lumpy","manta","delly"]

        for filename in filenames:

            try:
                # integrate with VCF version calculator
                annotate_vcf(filename, "GRCh37", "full")  # ANNOTATES VCF WITH ANNOTSV
            except Exception as e:
                logger.error("Error: {}".format(e))
                print("ANNOTATION ISSUE")
            try:
                data = read_json()  # READS BFF TEMPLATE
            except Exception as e:
                logger.error("Error: {}".format(e))
                print("READ BFF TEMPLATE ISSUE")
            try:
                pd = read_pandas(filename)  # USES PANDA TO READ ANNOTATED TSV
            except Exception as e:
                logger.error("Error: {}".format(e))
                print("READ ANNOTATED TSV ISSUE")
            try:
                parse_pandas(pd, data, filename)  # PARSES and WRITES ANNOTATION
            except Exception as e:
                logger.error("Error: {}".format(e))
                print("PARSE/ WRITE ISSUE")

            try:
                bff_post_processing(filename)
            except Exception as e:
                logger.error("Error: {}".format(e))
                print("POST PROCESSING ISSUE")


    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)
