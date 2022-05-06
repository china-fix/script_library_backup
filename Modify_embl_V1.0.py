'''
This script is design for modify the embl CDS features, and extract fixed length of the upstream seq
'''
# import glob
# import os
# import subprocess
import sys
import argparse
# import io
# from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq





def parse_args():
    parser=argparse.ArgumentParser(description="Welcome to use Xiao_Fei_Robot")
    parser.add_argument('--IN', required=True, type=str, metavar='FILENAME', help="the embl filename you want to modify")
    # parser.add_argument('--R_FOLDER', required=True, type=str, metavar='FILENAME', help="the complete reference genomes folder you want to blast to")
    parser.add_argument('--CUTLEN', default=100, type=float, metavar='DEFAULT 100', help="the cut length upstream the CDS")
    # parser.add_argument('--BLAST', default="blastn", type=str, metavar='command name', help="can be blastn(default) or tblastn")
    parser.add_argument('--OUT', default="out.embl", type=str, metavar='FILENAME', help="Output file name, out.embl as default")
    return parser.parse_args()

def main():
    args = parse_args()
    records = []
    features_up = []
    for seq_record in SeqIO.parse(args.IN, "embl"):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                if feature.location.strand == -1:
                    start_p = feature.location.end+0
                    end_p = feature.location.end + args.CUTLEN
                elif feature.location.strand == 1:
                    start_p = feature.location.start - args.CUTLEN
                    end_p = feature.location.start+0
                else:
                    start_p = feature.location.start - args.CUTLEN
                    end_p = feature.location.start+0
                feature_up_loc = FeatureLocation(start_p, end_p, strand= feature.location.strand)
                feature_up=SeqFeature(feature_up_loc, type='CDS')
                feature_up.qualifiers = feature.qualifiers
                try:
                    feature_up.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0]+'_UP'+str(args.CUTLEN)
                except:
                    pass
                features_up.append(feature_up)
        combine_features = seq_record.features + features_up
        seq_record.features = combine_features
        records.append(seq_record)
    SeqIO.write(records, args.OUT, "embl")







               









    
    
if __name__ == '__main__':
    sys.exit(main())
