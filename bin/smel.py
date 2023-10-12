#!/usr/bin/env python

import sys
import argparse

from .amplicon import download
from .amplicon import run

import warnings
warnings.filterwarnings(action='ignore')


def getOptions(argv):
    
    parser = argparse.ArgumentParser(
        prog='smel', 
        description='Automatic Microbiome Analysis Pipeline for SMEL')
    sub_parser = parser.add_subparsers(dest='command')
    
    # amplicon
    parser_amplicon = sub_parser.add_parser(
        'amplicon',
        help='Amplicon Analysis')
    sub_amplicon = parser_amplicon.add_subparsers(dest='subcommand')
    
    # download classifiers
    sub_amplicon_classifiers = sub_amplicon.add_parser('download',
                                                       help='Download amplicon classifiers')
    sub_amplicon_classifiers.add_argument(
        '-c', '--classifier',
        dest='classifier',
        metavar='STRING',
        type=str,
        help='Download amplicon classifiers (available: greengenes2, silva)'
    )
    
    # amplicon analysis options
    sub_amplicon = sub_amplicon.add_parser(
        'run', 
        help='Amplicon analysis with paired-end or PacBio CCS reads')
    sub_amplicon.add_argument(
        '-d', '--dirname',
        dest='dirname', 
        metavar='STRING [REQUIRED]', 
        type=str,
        help='Rawdata directory name to import')
    sub_amplicon.add_argument(
        '-m', '--metadata', 
        dest='metadata', 
        metavar='FILE [REQUIRED]', 
        type=str,
        help='Metadata file separated by tab')
    sub_amplicon.add_argument(
        '-c', '--classifier', 
        dest='classifier', 
        metavar='STR [OPTIONAL]', 
        type=str,
        help='Classifier path and name (greengenes2 or silva) [default : greengenes2]',
        default='greengenes2')
    sub_amplicon.add_argument(
        '-t', '--threads', 
        dest='threads', 
        metavar='INT [OPTIONAL]', 
        type=int,
        help='Number of threads during processing [default : 8]',
        default=8)
    sub_amplicon.add_argument(
        '--sampling_depth',
        dest='depth',
        metavar='INT [OPTIONAL]',
        type=int,
        help='Set initial sampling depth for rarefaction [default : 1000]',
        default=1000)
    sub_amplicon.add_argument(
        '--adapter_front',
        dest='adapt_fwd',
        metavar='SEQUENCE [OPTIONAL]',
        type=str,
        help="Forward adapter sequence (5'->3') [default: CCTACGGGNGGCWGCAG]",
        default='CCTACGGGNGGCWGCAG')
    sub_amplicon.add_argument(
        '--adapter_rev',
        dest='adapt_rev',
        metavar='SEQUENCE [OPTIONAL]',
        type=str,
        help="Reverse adapter sequence (5'->3') [default: GACTACHVGGGTATCTAATCC]",
        default='GACTACHVGGGTATCTAATCC')
    sub_amplicon.add_argument(
        '--ccs',
        dest='ccs',
        action='store_true',
        help='If PacBio CCS reads, add --ccs option')
    
    
    # shotgun
    parser_shotgun = sub_parser.add_parser(
        'shotgun',
        help='Shotgun metagenome analysis (Comming Soon)')
    sub_shotgun = parser_shotgun.add_subparsers(dest='subcommand')
    
    # paired end options
    sub_shotgun_paired = sub_shotgun.add_parser('paired', 
                                                help='Shotgun metagenome analysis with paired-end reads')
    sub_shotgun_paired.add_argument(
        '-d', '--dirname',
        dest='dirname', 
        metavar='DIRECTORY', 
        type=str,
        help='Rawdata directory name to import')
    
    
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    return args


def main():
    
    options = getOptions(sys.argv)
    command = options.command
    subcommand = options.subcommand
    
    # amplicon
    if command == 'amplicon':
        
        # download classifiers
        if subcommand == 'download':
            
            classifier = options.classifier
            download(classifier=classifier)
            
        
        # paired end
        elif subcommand == 'run':
            
            dirname = options.dirname
            metadata = options.metadata
            classifier = options.classifier
            depth = options.depth
            threads = options.threads
            adapter_front = options.adapt_fwd
            adapter_tail = options.adapt_rev
            ccs = options.ccs
            
            run(dirname=dirname,
                metadata=metadata,
                adapter_front=adapter_front,
                adapter_tail=adapter_tail,
                classifier=classifier,
                sampling_depth=depth,
                ccs=ccs,
                threads=threads)
        
    # shotgun
    elif command == 'shotgun':
        pass
    
    else:
        raise f'CommandNotFound: {subcommand}'

if __name__ == '__main__':
    main()