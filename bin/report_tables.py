#!/usr/bin/env python3

import re
import os
import sys
import base64
import shutil
import argparse
from glob import glob
from collections import defaultdict
from datetime import datetime

from lxml import html
from jinja2 import Environment, FileSystemLoader

# Parse plots to HTML file

date = datetime.strftime(datetime.now(), '%Y%m%d, %H:%M')

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--version', dest='version', help='WF version')
parser.add_argument('--doi', dest='doi', help='WF DOI') 
parser.add_argument('--templatedir', dest='templatedir', help='dir with template report')
args = parser.parse_args(sys.argv[1:])

# Get template, add bulma CSS stuff
templatefn = os.path.join(args.templatedir, 'report.html')
with open(os.path.join(args.templatedir, 'bulma.js')) as fp:
    bulma = fp.read()
with open(templatefn) as fp:
    template = Environment(loader=FileSystemLoader(args.templatedir)).from_string(fp.read())

# Plots from plotly
def get_plotly_html(fn):
    with open(fn) as fp:
        plothtml = html.parse(fp)
    plotbox = plothtml.find("body/div[@id='htmlwidget_container']")
    plotcode = plothtml.find("body/script[@type='application/json']")
    plot = html.tostring(plotbox, encoding=str).strip() + html.tostring(plotcode, encoding=str).strip().replace('\\n', '\\\\n').replace('\\\\\\n', '\\\\\\\\n')
    return plot

# Precursor plots
precplotfns = {'amount_precursors.html': '# of precursors',
        'missed_cleavages.html': '# of missed cleavages',
        }
precboxplots = {
        'retentiontime': 'Retention time',
        'precquant': 'Quantified precursors',
        'peakwidth': 'Peak width',
        'precerror': 'Precursor error',
        }

precplots = {}
plots = defaultdict(dict)

pdir = 'precplothtml'
nrchunks = len(glob(os.path.join(pdir, 'file__*__amount_precursors.html')))
for plotname in precplotfns:
    plot = defaultdict(list)
    for samfile in ['sample', 'file']:
        for chunk in range(nrchunks):
            pfile = os.path.join(pdir, f'{samfile}__{chunk}__{plotname}')
            if os.path.exists(pfile):
                plot[samfile].append({'title': precplotfns[plotname], 'plot': get_plotly_html(pfile)})
    precplots[plotname] = plot

for plotname in precboxplots:
    precplots[plotname] = []
    for chunk in range(nrchunks):
        plot = {}
        pfile = os.path.join(pdir, f'{plotname}__{chunk}.html')
        if os.path.exists(pfile):
            plot = {'title': precboxplots[plotname], 'plot': get_plotly_html(pfile)}
        precplots[plotname].append(plot)

# Protein/gene plots
featnames = [
       ('proteins', 'Proteins'), 
       ('genes', 'Genes'), 
]
featplotnames = [('nrfeats', 'Identifications'),
          ('missing_feats', 'Missing features'),
          ('ms1_quant', 'Quantities'),
          ('nn_quant', 'Quantities (no cross-run normalization)'),
          ('ms1nrprec', '# precursors with quant per protein'),
          ('ms1nrpep', '# peptides with MS2 quant per protein'),
          ]

featplotfns = {
        'missing_feats': ('missing_feats', False),
        'ms1_quant': ('quant', False),
        'nn_quant': ('nnquant', False),
        'ms1nrprec': ('nrp', False),
        # FIXME ALSO INCLUDE PEPTIDES
        'nrfeats': ('nrfeats', 'nrfeats__text.html'),
        }
featplots = defaultdict(dict)
for plotname, (pfn, textfn) in featplotfns.items():
    for featname, feattitle in featnames:
        featplots[plotname][featname] = []
        pdir = f'{featname}plots'
        for chunk in range(nrchunks):
            pfile = os.path.join(pdir, f'{pfn}__{chunk}.html')
            if os.path.exists(pfile):
                featplots[plotname][featname].append(get_plotly_html(pfile))
                if textfn:
                    with open(os.path.join(pdir, textfn)) as fp:
                        featplots[plotname][f'{featname}__text'] = fp.read().strip().split('\n')
    if all(x == [] for x in featplots[plotname].values()):
        featplots[plotname] = False

# JS libraries to add to HTML
# NB Only taking libs from one plot, same libs for other plots
for dirp, dirnames, fns in os.walk('.', followlinks=True):
    for fn in fns:
        if fn.endswith('.min.js'):
            src = os.path.join(dirp, fn)
            dst = os.path.join(dirp, fn.replace('.min.js', '.js'))
            shutil.copy(src, dst)
libs = []
for dirp, dirnames, fns in os.walk('precplothtml/sample__0__amount_precursors_files', followlinks=True):
    for fn in fns:
        srcfn = os.path.join(dirp, fn)
        if not os.path.exists(srcfn) or fn.endswith('.scss'):
            continue
        with open(srcfn) as fp:
            if fn.endswith('.js'):
                libs.append([fn, f'<script type="text/javascript">{fp.read()}</script>'])
            elif fn.endswith('.css') and not fn.endswith('.scss'):
                libs.append([fn, f'<style type="text/css">{fp.read()}</style>'])
if len(libs):
    with open('libs.js', 'w') as fp:
        for fn in ['htmlwidgets.js', 'plotly.js', 'typedarray.min.js', 'jquery.js', 'crosstalk.min.css',
                'crosstalk.js', 'plotly-htmlwidgets.css', 'plotly-latest.min.js']:
            try:
                lib = [x[1] for x in libs if x[0] == fn][0]
            except IndexError:
                print(f'Could not find library for HTML output {fn}, whic '
                        'should result from precursor table plotting')
                raise
            fp.write(f'{lib}\n')

# Summary table
# Titles for all tables
tabletitles = {
        'filename': 'Raw file', 
        'file': 'Raw file', 
        # TODO:
        #'no_pep_proteins': 'Peptides/protein (unique, median)',
        #'no_pep_genes': 'Peptides/gene (unique, median)',
        #'no_psm_proteins': 'PSMs/protein for quant (median)',
        #'no_psm_genes': 'PSMs/protein for quant (genecentric, median)',
        #'nr_proteins': '
        'Label': 'Sample',
        'nr_of_genes': 'Nr of genes ID',
        #'Non-shared (unique)': 'Peptides (unique, 1%FDR)',
        'precursorcount': '# precursors',
        'nr_scans': '# MS2 scans',
        }
summary_field_order = ['Label', 'precursorcount', 'nr_of_genes', 'nr_scans', 'perc_id',
        ]

summary_table = defaultdict(dict)
for fn in glob('counttable_qc.txt'):
    with open(fn) as fp:
        header = next(fp).strip('\n').split('\t')
        for line in fp:
            line = line.strip('\n').split('\t')
            fields = {f: line[ix] for ix, f in enumerate(header)}
            fname = fields.get('filename', fields.get('file'))
            for f in fields:
                summary_table[fname][f] = fields[f]
for _s, fields in summary_table.items():
    summary_fields = [x for x in summary_field_order if x in fields]
    break

# Missed cleavages table
miscleav = []
for chunk in range(nrchunks):
    with open(f'miscleav_qc.txt') as fp:
        head = next(fp).strip().split('\t')
        for line in fp:
            lnmap = {head[ix]: x for ix, x in enumerate(line.strip().split('\t'))}
            miscleav.append([lnmap['file'], lnmap['missed_cleavage'], lnmap['nrprec'], lnmap['percent']])

# Overlap table
overlap = defaultdict(dict)
for feattype, _ft in featnames:
    fn = f'{feattype}__overlap'
    if os.path.exists(fn):
        with open(fn) as fp:
            _h = next(fp)
            for line in fp:
                overlapnrfns, nr_feats = line.strip().split('\t')
                overlap[feattype][overlapnrfns] = nr_feats
    else:
        overlap[feattype] = False

## Warning box
#warnings = []
#for fn in glob('warnings*'):
#    if os.path.exists(fn):
#        with open(fn) as fp:
#            for line in fp:
#                if warn := line.strip():
#                    warnings.append(warn)
#

# Write to template
with open('report_groovy_template.html', 'w') as fp:
    fp.write(template.render(reportdate=date,
        reporttype='full',
        version=args.version,
        doi=args.doi,
        precplots=precplots,
        featplots=featplots,
        featnames=featnames,
        featplotnames=featplotnames,
        tabletitles=tabletitles,
        miscleav=miscleav,
        psmtables=False,
        summary_fields=summary_fields,
        summary_table=summary_table,
        overlap=overlap,
        #warnings=warnings,
        ))
with open('report_groovy_template_light.html', 'w') as fp:
    fp.write(template.render(reportdate=date,
        reporttype='light',
        version=args.version,
        doi=args.doi,
        precplots=precplots,
        featplots=featplots,
        featnames=featnames,
        featplotnames=featplotnames,
        tabletitles=tabletitles,
        miscleav=miscleav,
        psmtables=False,
        summary_fields=summary_fields,
        summary_table=summary_table,
        overlap=overlap,
        #warnings=warnings,
        ))
