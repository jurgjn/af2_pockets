
import collections
import difflib
import functools
import itertools
import math
import operator
import os
import random
import requests
import sys
import urllib
import urllib.parse
import yaml

import numpy as np
import scipy as sp
import scipy.stats
import scipy.special

import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

def uf(x):
    return '{:,}'.format(x)

def float_uk(s):
    return float(''.join(filter(lambda x: x in '0123456789.', s)))

def penc(s):
    """Encode exotic characters within identifiers using percent-encoding, e.g.:
        BF2.649 => BF2%2e649
        UH-232-(+) => UH-232-%28%2b%29
        bis(maltolato)oxovanadium(IV) => bis%28maltolato%29oxovanadium%28IV%29
    """
    #return urllib.parse.quote_plus(s) # Does not decode dots, e.g. BF2.649
    def penc_char(c):
        if c.isalnum() or c in ['_', '-']:
            return c
        else:
            return "%{0:0>2}".format(format(ord(c), "x"))

    return "".join(map(penc_char, s))

def pdec(s):
    """Opposite of pdec(). Sanity check:
    for id_ in ['BF2.649', 'UH-232-(+)', 'bis(maltolato)oxovanadium(IV)']:
        print(id_, '=>', penc(id_))
        assert id_ == pdec(penc(id_))
    """
    # unquote_plus seems to accept characters that aren't encoded by urllib.parse.quote_plus
    return urllib.parse.unquote_plus(s)

def pfile(asset_id=None, compound_id=None, struct_id=None, screen_id=None, step='{prev_steps}', base='{base}', suffix='', v=False):
    """
        Possible alternative to pf() above "asset id"
        For a small, finite set of "entity identifiers" (e.g. model; drug, model, cavity)
        define adhoc manual prefixes (using wildcard restrictions to ease parsing)
        - swiss_model_id:
        - swiss_protein_id:
    """
    if asset_id is None:
        l_asset_id = []
        if not(compound_id is None):
            if compound_id != '{}':
                (compound_ns, compound_name) = compound_id.split('_', maxsplit=1)
                if compound_ns == 'prism':
                    compound_pref = f'{compound_ns}_{compound_name[0]}'
                else:
                    compound_pref = compound_ns
            else:
                compound_pref = '{compound_pref}'
                compound_id = '{compound_id}'
            l_asset_id.append(compound_pref)
            l_asset_id.append(compound_id)

        if not(struct_id is None):
            if struct_id != '{}':
                struct_pref = os.path.join(struct_id[:2], struct_id[2:4], struct_id[4:6])
                #struct_pref = struct_id[:2]
            else:
                struct_pref = '{struct_pref}'
                struct_id = '{struct_id}'
            l_asset_id.append(struct_pref)
            l_asset_id.append(struct_id)

        if not(screen_id is None):
            l_asset_id.append(screen_id)

        asset_id = os.path.join(*l_asset_id)

    # add remaining elements: base, step, suffix
    filepath = os.path.join(base, step, f'{asset_id}{suffix}')
    if v: print('pfile: ', filepath)
    return filepath

def roc(labels, scores, *args, **kwargs):
    fpr, tpr, thresholds = sk.metrics.roc_curve(y_true=labels, y_score=scores, pos_label=True)
    plt.plot(fpr, tpr, *args, **kwargs)

def roc_ref(labels, scores, *args, **kwargs):
    fpr, tpr, thresholds = sk.metrics.roc_curve(y_true=labels, y_score=scores, pos_label=True)
    plt.plot(fpr, fpr, *args, **kwargs)

# CROC curves as defined in https://doi.org/10.1093/bioinformatics/btq140
def croc_transform(x, alpha=7):
    return (1 - math.exp(-alpha * x)) / (1 - math.exp(-alpha))

def croc(labels, scores, *args, **kwargs):
    fpr, tpr, thresholds = sk.metrics.roc_curve(y_true=labels, y_score=scores, pos_label=True)
    fpr_transform = list(map(croc_transform, fpr))
    plt.plot(fpr_transform, tpr, *args, **kwargs)

def croc_ref(labels, scores, *args, **kwargs):
    fpr, tpr, thresholds = sk.metrics.roc_curve(y_true=labels, y_score=scores, pos_label=True)
    fpr_transform = list(map(croc_transform, fpr))
    plt.plot(fpr_transform, fpr, *args, **kwargs)

def croc_xticks(xticks_=[0, 0.01, 0.05, 0.1, 0.2, 0.3, 1.0]):
    xticks_transform_ = list(map(croc_transform, xticks_))
    plt.gca().set_xticks(xticks_transform_)
    plt.gca().set_xticklabels(xticks_ )

def proc(labels = [True, False, False, True, False, True], scores = [1,2,3,1,2,3], *args, **kwargs):
    """ pseudo-ROC curves, e.g. Fig5 of Childs2019; Isik2015 Fig5
    proc([True, False, False, True], [1,2,3,4])
    proc([True, False, False, True, False, False, False, True], [1, 2, 3, 4, 1, 2, 3, 4])
    proc([True, True, True, True, False], [1, 2, 3, 4, 5])
    proc([False, False, True, False, False], [1, 2, 3, 4, 5])
    proc([False, True, False, True, False, True], [1, 2, 3, 4, 5, 6])
    proc([True, False, True, False, True, False], [1, 2, 3, 4, 5, 6])
    proc([True, False, True, False, True, False], [1, 1, 3, 4, 5, 6])

    Equal aspect ratio:
        plt.gca().set_aspect(1)
    """
    labels_scores_sorted = [(label, score) for label, score in sorted(zip(labels, scores), key=operator.itemgetter(1))]
    x = [0]
    y = [0]
    for label, score in itertools.groupby(labels_scores_sorted, key=operator.itemgetter(1)):
        x_next = x[-1] if len(x) > 0 else 0
        y_next = y[-1] if len(x) > 0 else 0
        for (label_, score_) in score:
            if label_:
                y_next += 1
            else:
                x_next += 1
        #x.append(x_next)
        x.append(score_)
        y.append(y_next)
    #plt.step(x, y, where='post') # Does not work in "tie cases" where equal scores have both labels
    plt.plot(x, y, *args, **kwargs)

def rrplot(labels = [True, False, False, True, False, True], scores = [1,2,3,1,2,3], *args, **kwargs):
    labels_scores_sorted = [(label, score) for label, score in sorted(zip(labels, scores), key=operator.itemgetter(1))]
    x = [0]
    y = [0]
    for label, score in itertools.groupby(labels_scores_sorted, key=operator.itemgetter(1)):
        x_next = x[-1] if len(x) > 0 else 0
        y_next = y[-1] if len(x) > 0 else 0
        for (label_, score_) in score:
            if label_:
                y_next += 1
            else:
                x_next += 1
        #x.append(x_next)
        x.append(score_)
        y.append(y_next)
    plt.step(x, y, where='pre', *args, **kwargs) # Does not work in "tie cases" where equal scores have both labels
