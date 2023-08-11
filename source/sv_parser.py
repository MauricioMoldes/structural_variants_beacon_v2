#!/usr/bin/env python

"""incomplete_submissions.py: Outputs a list of EGAS accessions and relevant information  for studies that have been published but not released."""

__author__ = "Arnau Soler, Mauricio Moldes"
__version__ = "0.1"
__maintainer__ = "Arnau Soler, Mauricio Moldes"
__email__ = "arnau.soler@crg.eu,mauricio.moldes@crg.eu"
__status__ = "test"

from pysam import VariantFile
import logging
import json

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
global position_endInteger
global position_startInteger
global position_refseqId
global position_assemblyId
global caseLevelData
global caseLevelData_zygosity
global caseLevelData_zigosity_label
global caseLevelData_zigosity_id
global caseLevelData_biosampleId
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


""" reads vcf """


def read_vcf():
    bcf_in = VariantFile("../data/delly.vcf.gz")  # auto-detect input format
    bcf_out = VariantFile('-', 'w', header=bcf_in.header)

    for rec in bcf_in.fetch('1', 10585, 14446418):
        bcf_out.write(rec)

""" write bff output in json array format """

def write_json():
    

if __name__ == "__main__":
    try:
        # configure logging
        log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s [in %(pathname)s:%(lineno)d]'
        logging.basicConfig(format=log_format)
        # execute main function
        read_vcf()
    except Exception as e:
        logger.error("Error: {}".format(e))
        sys.exit(-1)
