#!/usr/bin/env python3
"""Build two PDFs for the Swift-BAT E_peak (Task 1) and GBM-substitution (Task 2) work:
   - next_steps_figures.pdf      : figures + clear titles only
   - next_steps_explained.pdf    : same figures + titles + full explanations
"""
import os, matplotlib
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image, Table,
                                TableStyle, HRFlowable, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PImage
HERE=os.path.dirname(os.path.abspath(__file__))
fd=os.path.join(matplotlib.get_data_path(),"fonts/ttf")
pdfmetrics.registerFont(TTFont("DV",fd+"/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-Bold",fd+"/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-Obl",fd+"/DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFontFamily("DV",normal="DV",bold="DV-Bold",italic="DV-Obl",boldItalic="DV-Bold")
BODY=ParagraphStyle("b",fontName="DV",fontSize=10.5,leading=15,spaceAfter=7,alignment=TA_JUSTIFY)
H=ParagraphStyle("h",fontName="DV-Bold",fontSize=13,leading=16,spaceBefore=6,spaceAfter=5,textColor=colors.HexColor("#14365c"))
SUBH=ParagraphStyle("sh",fontName="DV-Bold",fontSize=10.5,leading=14,spaceBefore=6,spaceAfter=3,textColor=colors.HexColor("#333333"))
TITLE=ParagraphStyle("t",fontName="DV-Bold",fontSize=18,leading=22)
SUB=ParagraphStyle("s",fontName="DV-Obl",fontSize=11,textColor=colors.HexColor("#555555"))
CAP=ParagraphStyle("c",fontName="DV-Obl",fontSize=9.3,leading=12.5,textColor=colors.HexColor("#555555"),spaceBefore=3,spaceAfter=10)

def fig(path,w):
    iw,ih=PImage.open(path).size; im=Image(path,width=w*inch,height=w*inch*ih/iw); im.hAlign="CENTER"; return im
def box(txt):
    t=Table([[Paragraph(txt,BODY)]],colWidths=[6.7*inch])
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#eef3f9")),("BOX",(0,0),(-1,-1),0.6,colors.HexColor("#9ab6d6")),
        ("LEFTPADDING",(0,0),(-1,-1),10),("RIGHTPADDING",(0,0),(-1,-1),10),("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8)]))
    return t

FIG1="task1_bat_epeak.png"; FIG2="task2_gbm_substitution.png"
T1="Task 1 — Adding a Swift-measured prompt E_peak to the redshift estimator"
T2="Task 2 — Using Fermi GBM prompt variables in place of Swift's"
CAP1=("Left: the Swift BAT E_peak against redshift for the 255 sample GRBs with BAT data; the 57 directly-measured "
      "E_peak values (red) are shown separately from those derived from the photon index (light blue). The "
      "directly-measured E_peak is uncorrelated with redshift (r = +0.08, not significant). Right: the change in honest "
      "cross-validated r when BAT prompt variables are added to the estimator (baseline r = 0.47); error bars show the "
      "spread across cross-validation resamples. All changes are within the noise.")
CAP2=("The change in honest cross-validated r when the Swift peak flux is replaced by the Fermi GBM peak flux. On the "
      "original 65-GRB common sample it was +0.05 (within the noise); on a sample nearly doubled by a date+position "
      "crosswalk to the Fermi GBM catalog (116 GRBs, with a healthier baseline r = 0.41), it regresses to +0.01 and is "
      "positive in only 53% of resamples -- the small-sample hint does not hold up.")

# ---- explanation text (Task 1) ----
E1=[("What we ran and why",
     "Dr. Dainotti suggested obtaining an E_peak from the Swift BAT catalog itself -- from the 1-second peak spectrum, "
     "where the best-fit model is sometimes a cutoff power law -- so that we would not have to rely on Fermi. A key "
     "advantage is that this is a <b>measured Swift-side quantity for the whole sample</b>, so it avoids the circularity "
     "of the earlier Fermi transfer-learning approach (where an imputed E_peak was just a function of the Swift features "
     "the model already had)."),
    ("Data used",
     "We queried the HEASARC Swift BAT catalogs: <font face='DV-Obl'>swiftgrbba</font> (the BAT GRB catalog: 1-second "
     "peak photon index, 1-second peak flux, and the spectral-model flag) and <font face='DV-Obl'>swiftgrb</font> (which "
     "carries a direct <font face='DV-Obl'>bat_epeak</font>). Matching these to our 281-GRB redshift sample by GRB name "
     "-- including pre-2010 bursts that the BAT catalog names without the trailing letter (e.g. '050826' vs '050826A') "
     "-- gave <b>274 of 281 matched</b>, of which <b>255</b> have a 1-second peak photon index and peak flux and "
     "<b>57</b> have a directly-measured E_peak. Where no direct E_peak exists we estimated it from the BAT photon index "
     "via a relation calibrated on BAT bursts that have both: log10(E_peak) = -0.40*Gamma + 2.88 (correlation -0.45); "
     "this calibration never uses redshift, so the cross-validation stays honest."),
    ("What the results show",
     "The directly-measured BAT E_peak (57 GRBs) is <b>uncorrelated with redshift</b> (r = +0.08, not significant) -- "
     "confirming, now with a clean measured Swift quantity, that observer-frame E_peak carries no redshift signal "
     "(matching the Fermi result of r = 0.02). The BAT peak flux does correlate with redshift (r = -0.25) and the "
     "photon index weakly (r = -0.13), but adding any of the BAT prompt variables to the estimator (ridge on the "
     "winning-GLM feature set, baseline r = 0.47 on 255 GRBs) changes honest cross-validated r by only -0.0001 to "
     "-0.014 -- i.e. <b>no improvement</b>. The prompt peak flux's information is already captured by the existing "
     "Swift features. <b>Conclusion:</b> a directly-measured Swift prompt E_peak does not add redshift-predictive "
     "power, and neither do the other BAT prompt variables -- a clean, properly-powered, non-circular confirmation.")]
E2=[("What we ran and why",
     "The Swift peak flux is one of the seven input features of the redshift estimator. Dr. Dainotti's idea was that "
     "Fermi GBM measures the prompt emission more precisely than Swift's BAT, so feeding the model the GBM peak flux in "
     "place of the Swift one should give it a cleaner version of the <i>same</i> physical quantity and improve the "
     "redshift prediction. The experiment is therefore a <b>controlled head-to-head substitution</b>: we train the "
     "exact same model, on the exact same GRBs, twice -- once with the Swift peak-flux column, once with that single "
     "column swapped for the GBM peak flux -- and ask whether the cross-validated accuracy goes up. Nothing else about "
     "the model or the sample changes between the two runs, so any difference is attributable to the swap alone."),
    ("Data used, and how we enlarged the sample",
     "The experiment needs GRBs that have both a Fermi GBM peak flux and the Swift plateau features plus a redshift. "
     "Dr. Dainotti's hand-curated list had only <b>65</b> such bursts, which is too few to trust a +0.05. To get more, "
     "we matched the full 281-GRB redshift sample against the entire Fermi GBM catalog ourselves, using the fact that "
     "GRBs are named by date. For each of the 281 bursts we (1) read off its calendar day from its name (the YYMMDD "
     "stem), (2) pulled every Fermi GBM burst from that same day out of the HEASARC <font face='DV-Obl'>fermigbrst</font> "
     "catalog; (3) if exactly one Fermi burst occurred that day we matched it directly, and (4) if several did, we kept "
     "the one closest on the sky -- comparing the Swift position (from the <font face='DV-Obl'>swiftgrbba</font> catalog) "
     "with each Fermi position and requiring them to agree within 8 degrees, otherwise leaving the burst unmatched. This "
     "<b>date + sky-position crosswalk</b> raised the common sample to <b>122 GRBs</b> (116 with a usable GBM peak flux) "
     "-- nearly double the 65, and with a healthier baseline accuracy (r = 0.41 vs 0.20), so the enlarged result is the "
     "trustworthy one. The GBM peak flux is the Comptonized-fit photon flux (<font face='DV-Obl'>flnc_comp_phtflux</font>); "
     "both peak fluxes enter the model in log10."),
    ("The model and how we evaluated it",
     "The model is <b>ridge regression</b> on the seven Swift plateau features (log plateau flux and time, the plateau "
     "slope Alpha, column density N_H, photon index, peak flux and T90), each standardized to zero mean and unit "
     "variance first. Ridge is ordinary least-squares with an L2 penalty that shrinks the coefficients; its strength is "
     "chosen automatically by an internal cross-validation (RidgeCV). We use ridge rather than plain least-squares "
     "because these common samples are small (65-116 GRBs) and the features are correlated, which makes an unpenalized "
     "fit unstable -- ridge keeps the comparison stable and fair. We score it with <b>honest repeated 5-fold "
     "cross-validation</b>: the GRBs are split into 5 folds, the model is trained on 4 and predicts the held-out fifth "
     "(so every burst is predicted by a model that never saw it), and this is repeated 400 times with different random "
     "splits. For each split we take the Pearson r between predicted and true redshift (converted back to linear z), "
     "<b>paired</b> so that the Swift and GBM models see identical training points and folds; the reported change is the "
     "average of (GBM r &minus; Swift r) across resamples, and '% positive' is the fraction of resamples in which the "
     "GBM version actually won."),
    ("What the results show",
     "On the original 65 GRBs, replacing the Swift peak flux with the GBM peak flux raised r by about +0.05 (0.20 -> "
     "0.25) and the GBM peak flux correlated slightly more with redshift (-0.34 vs -0.29) -- suggestive, consistent with "
     "Dr. Dainotti's expectation that the GBM measurement is more precise, but within the noise. <b>On the enlarged "
     "116-GRB sample</b> (healthier baseline r = 0.41), the same substitution gives only <b>+0.01 ± 0.09, positive in "
     "53% of resamples -- essentially a coin flip</b>, and the GBM peak flux's correlation advantage shrinks. "
     "<b>Conclusion:</b> the small-sample +0.05 was noise; with roughly double the data the GBM peak flux does not "
     "improve redshift prediction.")]

def build(path, explained):
    doc=SimpleDocTemplate(path,pagesize=letter,leftMargin=0.9*inch,rightMargin=0.9*inch,topMargin=0.8*inch,bottomMargin=0.75*inch)
    S=[Paragraph("GRB Redshift Estimator: Next-Steps Results",TITLE),Spacer(1,4),
       Paragraph(("Swift BAT E_peak and Fermi GBM substitution "+("(with explanations)" if explained else "(figures)"))+" · October 2026",SUB),
       Spacer(1,8),HRFlowable(width="100%",color=colors.HexColor("#cccccc")),Spacer(1,8)]
    S+=[box("<b>Summary.</b> Two follow-up directions from Dr. Dainotti. <b>Task 1</b> gets a measured E_peak from the "
            "Swift BAT catalog and adds it to the estimator (no improvement -- a clean, non-circular confirmation that "
            "prompt E_peak carries little redshift signal). <b>Task 2</b> substitutes the Fermi GBM peak flux for the "
            "Swift one: a small-sample +0.05 did not survive a crosswalk-enlarged sample (it regressed to +0.01), so it "
            "is not a real gain. Neither task improves the estimator."),Spacer(1,12)]
    for T,FIGP,CAPT,EXPL,w in [(T1,FIG1,CAP1,E1,6.6),(T2,FIG2,CAP2,E2,6.0)]:
        block=[Paragraph(T,H),Spacer(1,3),fig(FIGP,w),Paragraph("<i>"+CAPT+"</i>",CAP)]
        S+=[KeepTogether(block)]
        if explained:
            for sh,txt in EXPL: S+=[Paragraph(sh,SUBH),Paragraph(txt,BODY)]
        S+=[Spacer(1,14)]
    doc.build(S); print("built",os.path.basename(path))

build(HERE+"/next_steps_figures.pdf", explained=False)
build(HERE+"/next_steps_explained.pdf", explained=True)
