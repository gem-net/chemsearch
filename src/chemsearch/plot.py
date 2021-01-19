"""Plot rdkit molecule object."""
import os
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import transforms, text
from rdkit.Chem import Draw


FIGSIZE = 200
MAX_FONTSIZE = 12
_logger = logging.getLogger(__name__)


def plot_mol(mol):
    fig = Draw.MolToMPL(mol, size=(FIGSIZE, FIGSIZE))
    ax = fig.get_axes()[0]
    #     ax.set_position([0, 0, 1, 1])
    ax_to_disp = ax.transAxes
    disp_to_fig = fig.transFigure.inverted()
    bb_dl = ax.dataLim
    p0, p1 = bb_dl.p0, bb_dl.p1
    n0, n1 = disp_to_fig.transform(ax_to_disp.transform([p0, p1]))
    bb = mpl.transforms.Bbox([n0, n1])
    l, b, r, t = bb_dl.extents
    ax.set_xlim(l, r)
    ax.set_ylim(b, t)
    ax.set_position(bb)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.set_frame_on(False)
    annots = [i for i in ax.get_children() if type(i) == mpl.text.Annotation]
    for annot in annots:
        if annot.get_fontsize() > MAX_FONTSIZE:
            annot.set_fontsize(MAX_FONTSIZE)
    return fig, ax


def save_svg_if_not_present(mol, svg_path):
    if os.path.exists(svg_path):
        return
    hf, ax = plot_mol(mol)
    hf.savefig(svg_path, bbox_inches='tight')
    plt.close(hf)
