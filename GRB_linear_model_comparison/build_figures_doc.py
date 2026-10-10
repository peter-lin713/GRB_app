#!/usr/bin/env python3
"""Compile the key audit figures into one PDF."""
import os, matplotlib
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image, Table,
                                TableStyle, HRFlowable, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PImage

fd = os.path.join(matplotlib.get_data_path(), "fonts/ttf")
pdfmetrics.registerFont(TTFont("DV", fd+"/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-Bold", fd+"/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-Obl", fd+"/DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-Bold", italic="DV-Obl", boldItalic="DV-Bold")

doc = SimpleDocTemplate("redshift_estimator_figures.pdf", pagesize=letter,
                        leftMargin=0.9*inch, rightMargin=0.9*inch, topMargin=0.8*inch, bottomMargin=0.75*inch)
BODY = ParagraphStyle("b", fontName="DV", fontSize=10.5, leading=15, spaceAfter=6, alignment=TA_JUSTIFY)
H = ParagraphStyle("h", fontName="DV-Bold", fontSize=12, leading=15, spaceBefore=4, spaceAfter=5,
                   textColor=colors.HexColor("#14365c"))
CAP = ParagraphStyle("c", fontName="DV", fontSize=9.3, leading=12.5, spaceBefore=4, spaceAfter=2,
                     textColor=colors.HexColor("#333333"), alignment=TA_JUSTIFY)
TITLE = ParagraphStyle("t", fontName="DV-Bold", fontSize=19, leading=23)
SUB = ParagraphStyle("s", fontName="DV-Obl", fontSize=11, textColor=colors.HexColor("#555555"))

def fig(path, w_in):
    iw, ih = PImage.open(path).size
    im = Image(path, width=w_in*inch, height=w_in*inch*ih/iw); im.hAlign = "CENTER"; return im

FIGS = [
 ("comparison_linear_vs_glm.png", 4.7,
  "1.  Winning GLM and a plain linear model, under the flagship protocol and under honest cross-validation",
  "Predicted vs. observed redshift. Columns: the winning GLM vs. a plain 7-feature linear model. "
  "Top row: the flagship protocol (M-estimator outlier removal + cross-validation + 2&#963; clip). "
  "Bottom row: honest cross-validation with nothing removed. Under the flagship protocol the two models "
  "give similar results (r &#8776; 0.63&#8211;0.66); under honest cross-validation both are near r &#8776; "
  "0.46&#8211;0.53. The difference between the two sets of numbers tracks the evaluation protocol rather than "
  "the choice of model."),
 ("true_vs_pred_sidebyside.png", 6.7,
  "2.  Two ways of removing 5%: the M-estimator (residual-based) and a luminosity cut",
  "Both panels use the identical honest out-of-fold GLM predictions; only the highlighted 5% differs. "
  "Left: the M-estimator removes points that lie off the 1:1 line (the worst-predicted), and r changes "
  "0.53 &#8594; 0.62. Right: a luminosity cut removes the faintest 5%, which sit spread around the line, and "
  "r changes little (0.53 &#8594; 0.54). The increase in r is associated specifically with removing the "
  "worst-predicted points."),
 ("outlier_train_vs_eval.png", 6.7,
  "3.  Training vs. evaluation contributions of the removed points",
  "The 15 M-estimator-removed points decomposed into a training effect vs. an evaluation effect. Removing them "
  "from the test set changes r by +0.088; removing them from training only changes it by +0.021. Right panel: "
  "a model trained on the 266 remaining bursts predicts the 15 held-out points with about 2.8&#215; the error "
  "(r &#8776; &#8722;0.22) &mdash; they are difficult to predict, so most of the effect comes from whether or "
  "not they are in the test set."),
 ("luminosity_cut_train_vs_eval.png", 4.3,
  "4.  The same decomposition for a luminosity-based (physical) 5% cut",
  "The same train/test decomposition, but removing the faintest 5% by luminosity (a Malmquist-style cut) "
  "instead of the M-estimator&#8217;s residual-based 5%. Both the training and evaluation effects are near "
  "zero (&#177;0.003); the change in r seen in Figure 3 is associated with the residual-based selection rather "
  "than with removing faint bursts."),
 ("heldout_r_distribution.png", 5.9,
  "5.  Held-out performance with the M-estimator applied to the training set only",
  "Distribution of honest held-out r over 1000 random 20-point holdouts, with the M-estimator fit on the "
  "training split only and the GLM retrained on it. The mean (0.52), median (0.54) and mode (0.64) are marked, "
  "along with the value obtained when points are removed from the full sample (0.63). The left-skew explains "
  "why a typical holdout lands near 0.6 while the mean is ~0.52; the mode is close to the full-sample value."),
 ("outlier_reversed_relation.png", 5.0,
  "6.  Flux&#8211;redshift relation of the removed bursts compared with the rest",
  "Plateau flux vs. redshift. The bulk (blue) shows the usual weak anti-correlation (fainter at higher z); "
  "the removed points (red) show the opposite trend &mdash; brighter at higher z (r = +0.69). A single global "
  "model learns the bulk&#8217;s slope and so predicts these points in the opposite direction, which is why the "
  "M-estimator identifies them as outliers."),
 ("luminosity_vs_redshift.png", 5.0,
  "7.  The removed bursts in the luminosity&#8211;redshift plane",
  "Plateau luminosity vs. redshift with the flux-limit curve. The M-estimator-removed bursts (red) fall at the "
  "extremes of the L&#8211;z sequence &mdash; faint-and-nearby or bright-and-distant &mdash; rather than along a "
  "single feature axis. These are among the objects a single relation predicts least well."),
 ("heldout_xray_optical_remove10.png", 5.3,
  "8.  Held-out performance on the X-ray + optical sample, with and without the difficult points",
  "Distribution of honest held-out r over 1000 random 80/20 splits on the combined X-ray + optical sample "
  "(296 bursts), with the M-estimator applied to the training split only and the GLM retrained on it. "
  "Top: all bursts (mean 0.48, median 0.49, mode 0.51). Bottom: with the 10 difficult-to-predict points removed "
  "beforehand (286 bursts; mean 0.54, median 0.54, mode 0.54). Removing them shifts the distribution up by about "
  "0.06 and removes the left-skew seen in the top panel. The red dashed line marks the flagship value (difficult "
  "points removed from the full sample, then cross-validated): r = 0.56."),
 ("heldout_xray_only_remove10.png", 5.3,
  "9.  The same procedure, restricted to the X-ray-only sample",
  "As in Figure 8, restricted to the X-ray-only subset (223 bursts). Of the 10 difficult points, 8 are X-ray; "
  "the other two (100508A and 160314A) are optical and are not in this subset. Top: all 223 (mean 0.51, median "
  "0.52, mode 0.54). Bottom: with the 8 X-ray difficult points removed (215 bursts; mean 0.57, median 0.57, mode "
  "0.56). The shift is again about +0.06. The red dashed line marks the flagship value (points removed from the "
  "full sample, then cross-validated): r = 0.59."),
]

S = [Paragraph("GRB Redshift Estimator: Analysis Figures", TITLE), Spacer(1,4),
     Paragraph("Figures from the analysis of the redshift estimator &#183; October 2026", SUB),
     Spacer(1,8), HRFlowable(width="100%", color=colors.HexColor("#cccccc")), Spacer(1,6)]
box = Table([[Paragraph("<b>About these figures.</b> These figures examine the reported r &#8776; 0.71 and how it "
    "depends on the evaluation protocol. Under an honest cross-validation with no point removal, a plain linear "
    "model and the winning GLM both give r &#8776; 0.5. The figures compare the two models, the two kinds of 5% "
    "cut (residual-based vs. luminosity-based), the training-vs-evaluation decomposition, the honest held-out "
    "distribution, and the properties of the removed bursts.", BODY)]],
    colWidths=[6.7*inch])
box.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#eef3f9")),
    ("BOX",(0,0),(-1,-1),0.6,colors.HexColor("#9ab6d6")),("LEFTPADDING",(0,0),(-1,-1),10),
    ("RIGHTPADDING",(0,0),(-1,-1),10),("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8)]))
S += [box, Spacer(1,10)]

for path, w, head, cap in FIGS:
    S.append(KeepTogether([Paragraph(head, H), Spacer(1,3), fig(path, w),
                           Paragraph("<i>"+cap+"</i>", CAP), Spacer(1,14)]))

doc.build(S)
print("built redshift_estimator_figures.pdf")
