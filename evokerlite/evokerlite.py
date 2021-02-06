#!/usr/bin/env python
from .batches import Batches
from .snp_posterior import SNPPosterior
from .variants import Variants
from .samples import Samples
from .genotypes import Genotypes
from .intensity import TextIntensity, BinaryIntensity
import logging
from argparse import ArgumentParser
from collections import defaultdict
import shutil
import os
import re
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
logging.basicConfig(level=logging.INFO)


ALLOWED_IMAGE_FORMATS = ('.pdf', '.png')
DEFAULT_IMAGE_FORMAT = '.png'

# UKB chromosome convention: X=23,Y=24,XY=25,MT=26
UKB_SEX_CHROMOSOMES = {
    '23': 'X',
    '24': 'Y',
    '25': 'XY',
    'X': 'X',
    'Y': 'Y',
    'XY': 'XY',
}

RE_FAM = re.compile(r'.*\.fam$')
RE_BIM = re.compile(r'.*\.bim$')
RE_UKB_BIM = re.compile(r'ukb_snp_chr([0-9,X,x,Y,y]+)_v2\.bim$')


def file_format_check(ext, fallback=None):
    if ext not in ALLOWED_IMAGE_FORMATS:
        raise Exception('File format must be on one of either: {}'.format(
            ', '.join(ALLOWED_IMAGE_FORMATS)))


def scatter(ax, t, m):
    ax.scatter(t[:, 0], t[:, 1],
               color=m['color'],
               s=50,
               alpha=m.get('alpha', 0.5),
               lw=m.get('lw', 0.25),
               edgecolor='black',
               label=m['label']
               )


def setup_ax(ax, title, transform):
    ax.set_title(title, fontsize=14, weight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if transform:
        xlabel = 'Contrast: log2(A/B)'
        ylabel = 'Strength: log2(A*B)/2'
    else:
        xlabel = 'A'
        ylabel = 'B'
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()


class EvokerLite:

    def __init__(self, bfile_path=None, int_path=None, fam_path=None,
                 bnt_path=None, bim_path=None, bed_path=None,
                 snp_posterior_batch_path=None, snp_posterior_path=None,
                 exclude_list_path=None, ukbiobank=False):

        self.ukbiobank = ukbiobank
        self.samples = Samples(fam_path or (bfile_path + '.fam'), ukbiobank)
        self.variants = Variants(bim_path or (bfile_path + '.bim'))
        n_samples = self.samples.get_n_samples()
        n_variants = self.variants.get_n_variants()
        self.genotypes = Genotypes(bed_path or (bfile_path + '.bed'),
                                   n_samples, n_variants)

        if ukbiobank:
            if snp_posterior_path:
                if not snp_posterior_batch_path:
                    raise Exception('When plotting snp posterior\'s in UK'
                                    'Biobank requires a batch file.')
                self.snp_posterior_batches = Batches(snp_posterior_batch_path)
                logging.debug(snp_posterior_path)
                self.snp_posterior = SNPPosterior(file_path=snp_posterior_path,
                                                  n_batches=self.snp_posterior_batches.get_n_batches())
        if bnt_path:
            self.intensities = BinaryIntensity(
                bnt_path, self.samples, ukbiobank)
        elif int_path:
            self.intensities = TextIntensity(
                int_path, self.variants, self.samples)
        else:
            bnt_path = bfile_path + '.bnt'
            int_path = bfile_path + '.int'
            if os.path.isfile(bnt_path):
                self.intensities = BinaryIntensity(
                    bnt_path, self.samples, ukbiobank)
            elif os.path.isfile(int_path):
                self.intensities = TextIntensity(
                    int_path, self.samples, ukbiobank)
            else:
                raise Exception('Must provide a binary or intensity file')

        if exclude_list_path:
            exclude_list = []
            with open(exclude_list_path) as f:
                for line in f:
                    sample = line.strip()
                    exclude_list.append(sample)
            self.exclude_list = exclude_list

    def plot(self, variant_name, batch=None, ellipses=None, ax=None, transform=None):
        variant_index = self.variants.get_index(variant_name)
        variant = self.variants.get_variant(variant_index)
        A1 = variant.get_A1()
        A2 = variant.get_A2()
        chrom = variant.get_chrom()

        genotype_mapping = (
            {'code': '00', 'name': 'homozygous A1', 'color': 'blue',
                'label': '{A1}{A1}'.format(A1=A1, A2=A2)},
            {'code': '10', 'name': 'heterozygous', 'color': 'limegreen',
                'label': '{A1}{A2}|{A2}{A1}'.format(A1=A1, A2=A2)},
            {'code': '11', 'name': 'homozygous A2', 'color': 'red',
                'label': '{A2}{A2}'.format(A1=A1, A2=A2)},
            {'code': '01', 'name': 'missing', 'color': 'grey',
                'alpha': 0.9, 'lw': 0.75, 'label': 'No call'},
        )
        genotypes = self.genotypes.get_genotypes(variant_index)
        if self.ukbiobank and transform is not False:
            transform = True
        xy = self.intensities.get_intensities_for_variant(
            variant_index, transform)

        logging.debug(xy)
        plot_sexes = False
        if self.ukbiobank:
            if batch is None:
                raise Exception(
                    'Must specify a batch when plotting UK Biobank data.')
            sample_batches = self.samples.get_batches()
            batch_indices = sample_batches == batch
            logging.debug(batch_indices.shape)
            logging.debug(genotypes.shape)
            genotypes = genotypes[batch_indices]
            xy = xy[batch_indices]

            if chrom in UKB_SEX_CHROMOSOMES.keys():
                plot_sexes = True
                sexes = self.samples.get_sex()
                sexes = sexes[batch_indices]
                male_indices = sexes == 1
                female_indices = sexes == 2
                ellipses = None
            if ellipses:
                if plot_sexes:
                    batch_index = self.snp_posterior_batches.get_index(batch)
                    ellipses = self.snp_posterior.get_batch_ellipse_points(
                        variant_index, batch_index)
                else:
                    ellipses = None
        if ax:
            fig = False
        else:
            if plot_sexes:
                fig, ax = plt.subplots(1, 2, figsize=(20, 10))
            else:
                fig, ax = plt.subplots(figsize=(10, 10))
        for m in genotype_mapping:
            code = m['code']
            if plot_sexes:
                t = xy[(genotypes == code) & male_indices]
                scatter(ax[0], t, m)
                t = xy[(genotypes == code) & female_indices]
                scatter(ax[1], t, m)
            else:
                t = xy[genotypes == code]
                scatter(ax, t, m)
            if ellipses and code in ellipses:
                ellipse_points = ellipses[code]
                if plot_sexes:
                    ax[0].plot(ellipse_points[:, 0],
                               ellipse_points[:, 1], color='black')
                    ax[1].plot(ellipse_points[:, 0],
                               ellipse_points[:, 1], color='black')
                else:
                    ax.plot(ellipse_points[:, 0],
                            ellipse_points[:, 1], color='black')

        title = variant_name
        if self.ukbiobank:
            c = UKB_SEX_CHROMOSOMES.get(chrom, chrom)
        else:
            c = chrom
        title += ' | chr{}'.format(c)
        if batch:
            title += ' | {}'.format(batch)
        title += ' | n={}'.format(len(genotypes))
        if plot_sexes:
            setup_ax(ax[0], title + ' | Male', transform)
            setup_ax(ax[1], title + ' | Female', transform)
        else:
            setup_ax(ax, title, transform)

        if fig:
            return fig

    def save_plot(self, variant_name, outpath=None, batch=None, transform=None, ellipses=None, image_format=DEFAULT_IMAGE_FORMAT):
        fig = self.plot(variant_name=variant_name, batch=batch,
                        ellipses=ellipses, transform=transform)
        image_format = image_format.lower()
        file_format_check(image_format)
        if batch:
            filename = '{}_{}{}'.format(variant_name, batch, image_format)
        else:
            filename = variant_name + image_format
        if not outpath:
            outpath = filename
        elif os.path.isdir(outpath):
            outpath = os.path.join(outpath, filename)
        else:  # The only other option is a filename
            filename, ext = os.path.splitext(outpath)
            ext = ext.lower()
            if not ext:
                outpath = filename + image_format
            else:
                file_format_check(ext)
                outpath = filename
        print('Saving to: {}'.format(outpath))
        fig.savefig(outpath)
        plt.close('all')

    def save_all_batches(self, variant_name, outdirectory, transform=None, ellipses=None, image_format=DEFAULT_IMAGE_FORMAT):
        if not os.path.exists(outdirectory):
            os.mkdir(outdirectory)
        directory = os.path.join(outdirectory, variant_name)
        if os.path.exists(directory):
            print('Warning: the directory to hold all of the plots ({}) already exists.'.format(
                directory))
        else:
            os.mkdir(directory)
        batches = sorted(list(set(self.samples.get_batches())))
        for batch in batches:
            self.save_plot(variant_name, outpath=directory,
                           batch=batch, ellipses=ellipses, transform=transform)


def get_rsids(bim):
    rsids = []
    print(bim)
    with open(bim) as f:
        for line in f:
            tokens = line.split()
            rsids.append(tokens[1])
    return set(rsids)


def get_chrom_to_rsids(rsids, data):
    rsids = set(rsids)
    ls = os.listdir(data)
    def dd(x): return os.path.join(data, x)
    chrom2rsids = defaultdict(set)
    found = set()
    for filename in ls:
        m = RE_UKB_BIM.match(filename)
        if m:
            chrom = m.groups()[0]
            bim_rsids = get_rsids(dd(filename))
            overlap = bim_rsids & rsids
            if overlap:
                logging.debug(chrom)
                logging.debug(overlap)
                chrom2rsids[chrom].update(overlap)
                found.update(overlap)
                if found == rsids:
                    break
    return chrom2rsids


def plot_uk_biobank(data, output, rsids, transform=True, snp_posterior=False, fam=None):
    chrom2rsids = get_chrom_to_rsids(data, rsids)
    ls = os.listdir(data)

    if fam:
        famfile = fam
    else:
        for filename in ls:
            m = RE_FAM.match(filename)
            if m:
                famfile = filename
                break
        else:
            raise Exception('Directory should contain a single fam file.')

    def dd(x): return os.path.join(data, x)

    for chrom, _rsids in chrom2rsids.items():
        params = {
            'bed_path': dd('ukb_cal_chr{}_v2.bed'.format(chrom)),
            'fam_path': dd(famfile),
            'bim_path': dd('ukb_snp_chr{}_v2.bim'.format(chrom)),
            'bnt_path': dd('ukb_int_chr{}_v2.bin'.format(chrom)),
            'ukbiobank': True,
        }
        if snp_posterior:
            params['snp_posterior_path'] = dd(
                'ukb_snp_posterior_chr{}.bin'.format(chrom))
            params['batch_path'] = dd('ukb_snp_posterior.batch')
        el = EvokerLite(**params)
        for rsid in _rsids:
            el.save_all_batches(variant_name=rsid, outdirectory=output,
                                transform=transform, ellipses=snp_posterior)


def get_file_prefix_to_rsids(rsids, data):
    rsids = set(rsids)
    def dd(x): return os.path.join(data, x)
    file_prefix2rsids = defaultdict(set)
    found = set()
    for filename in os.listdir(data):
        m = RE_BIM.match(filename)
        if m:
            file_prefix = os.path.splitext(filename)[0]
            bim_rsids = get_rsids(dd(filename))
            overlap = bim_rsids & rsids
            if overlap:
                logging.debug(file_prefix)
                logging.debug(overlap)
                file_prefix2rsids[file_prefix].update(overlap)
                found.update(overlap)
                if found == rsids:
                    return file_prefix2rsids
    else:
        raise Exception(
            f'Only found {len(found)}/{len(rsids)} rsids in data. Aborting.')


def plot(data, output, rsids, fam=None):
    file_prefix_to_rsids = get_file_prefix_to_rsids(rsids, data)
    def dd(x): return os.path.join(data, x)

    params = {}
    if fam:
        params['fam_path'] = dd(fam)
    for file_prefix, _rsids in file_prefix_to_rsids.items():
        params = {
            'bfile_path': dd(file_prefix)
        }
        el = EvokerLite(**params)
        for rsid in _rsids:
            el.save_plot(variant_name=rsid, outpath=output)


def cli():
    parser = ArgumentParser()
    parser.add_argument('--ukb',
                        action='store_true',
                        help='flag to indicate data is in UK Biobank format'
                        )
    parser.add_argument('-d', '--data',
                        type=str,
                        required=True,
                        help='directory of PLINK/intensity data'
                        )
    parser.add_argument('-f', '--fam',
                        type=str,
                        help='location of fam file (specify if not in the same directory as of PLINK/instensity data',
                        )
    parser.add_argument('-o', '--output',
                        type=str,
                        help='directory to save plots',
                        )
    parser.add_argument('-r', '--rsids',
                        type=str,
                        required=True,
                        help='text file with rsids to create cluster plots from',
                        )
    parser.add_argument('--no-transform',
                        action='store_true',
                        help='flag to not plot UKBiobank data in contrast/strength coordinates'
                        )
    parser.add_argument('--snp-posterior',
                        action='store_true',
                        help='plot UKBiobank SNP Posterior'
                        )

    args = parser.parse_args()

    output = args.output
    if not output:
        output = os.getcwd()

    with open(args.rsids) as f:
        rsids = [x.strip() for x in f.readlines()]

    if args.ukb:
        transform = not args.no_transform
        if not transform:
            snp_posterior = False
        else:
            snp_posterior = args.snp_posterior
        plot_uk_biobank(args.data, output, rsids,
                        transform, snp_posterior, args.fam)
    else:
        plot(args.data, output, rsids, args.fam)


if __name__ == '__main__':
    cli()
