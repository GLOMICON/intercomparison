#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Modifs : - passage en python 3 (Greg and Loh; 12/11/2019)
#          - Update du code, ajout de compteurs, argumentParser instead of optionParser (Greg; 16/10/2020)
#          - Fixed indentation bug, Reorganised arguments, Default behavior is not whole-word-only substitution.
#            Added a 'pattern_search' argument to substitute anything (Greg; 27/04/2021)

import logging
import re
import sys
from argparse import ArgumentParser


def parse_substitution_table(arguments):
    logging.info(f'Step 1. Reading substitution table ({arguments.substitution_table})')
    sub_table = {}
    sub_counts = {}

    with open(arguments.substitution_table) as sub:
        for line in sub:
            line = line.strip()
            parsed = line.split('\t')
            if parsed[0] not in sub_table and len(parsed) == 2:
                if arguments.reverse:
                    sub_table[parsed[1]] = parsed[0]
                    sub_counts[parsed[1]] = 0
                else:
                    sub_table[parsed[0]] = parsed[1]
                    sub_counts[parsed[0]] = 0

            else:
                logging.error(f'Following substitution pattern is not recognised (not 2 columns) or duplicated:\n'
                              f'    {parsed[0]}')

    return sub_table, sub_counts


def process_substitution(arguments, sub_table, sub_counts):
    logging.info(f'Step 2. Performing substitution (outfile={arguments.output_file})')
    if arguments.pattern_search:
        logging.warning('Option --pattern_search is enabled, substitutions will happen within words.')
    else:
        logging.info('Whole-word substitutions only. Make sure that it\'s what you need.')

    with open(arguments.input_file) as file, open(arguments.output_file, 'w') as out:
        for line in file:
            line = line.strip()
            for motif, subst in sub_table.items():
                if motif in line:
                    sub_counts[motif] += 1
                    if arguments.pattern_search:
                        line = re.sub(motif.replace('|', '\\|'), subst, line)
                    else:
                        line = re.sub(rf'\b{motif}\b', subst, line)
                        # \b stands for 'word boundary'

            out.write(f'{line}\n')

    logging.info(f'Number of patterns to substitute: {len(sub_counts)}')
    once_or_more = [k for k in sub_counts if sub_counts[k] > 0]
    logging.info(f'Number of patterns found at least once: {len(once_or_more)}')
    twice_or_more = sorted([(k, sub_counts[k]) for k in sub_counts if sub_counts[k] > 1],
                           key=lambda _: _[1], reverse=True)
    if len(twice_or_more) > 0:
        pack = '\n'.join([f'    {k} ({ct})' for k, ct in twice_or_more])
        logging.warning(f'Patterns substituted more than once:\n{pack}')


if __name__ == '__main__':
    # # OPTION PARSER # #
    parser = ArgumentParser(usage='name_substituer.py -i input_file -o output_file -s substitution_table (-v)')
    parser.add_argument('-i', '--input', dest='input_file', help='Name for input file')
    parser.add_argument('-s', '--subtable', dest='substitution_table',
                        help='Name for substitute 2-columns tabulated file [motif_to_replace\treplacement_motif]')
    parser.add_argument('-o', '--output', dest='output_file', help='Name for output file', required=True)
    parser.add_argument('-p', '--pattern', dest='pattern_search', help='Search for pattern instead of full word',
                        action='store_true', default=False)
    parser.add_argument('-r', '--reverse', dest='reverse', help='Reverse substitution', action='store_true',
                        default=False)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if args.input_file is None or args.substitution_table is None:
        logging.error('Missing critical arguments -i and/or -s. For help, use -h. Exiting.')
        sys.exit(-1)

    logging.info(f'Running name_substituer.py on {args.input_file}.')
    substitution_table, substitution_counts = parse_substitution_table(args)
    process_substitution(args, substitution_table, substitution_counts)

