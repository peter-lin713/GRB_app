#!/usr/bin/env python3
"""Build the transfer-learning feasibility note as a PDF."""
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, ListFlowable, ListItem, HRFlowable)

import os, matplotlib
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
_fd = os.path.join(matplotlib.get_data_path(), "fonts/ttf")
pdfmetrics.registerFont(TTFont("DV", _fd+"/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-Bold", _fd+"/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-Obl", _fd+"/DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFont(TTFont("DVMono", _fd+"/DejaVuSansMono.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-Bold", italic="DV-Obl", boldItalic="DV-Bold")

doc = SimpleDocTemplate("transfer_learning_feasibility.pdf", pagesize=letter,
                        leftMargin=0.9*inch, rightMargin=0.9*inch,
                        topMargin=0.8*inch, bottomMargin=0.8*inch)
ss = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=ss["Normal"], fontName="DV", fontSize=10.5,
                      leading=15, spaceAfter=7, alignment=TA_JUSTIFY)
H1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName="DV-Bold", fontSize=13,
                    leading=16, spaceBefore=15, spaceAfter=5, textColor=colors.HexColor("#14365c"))
TITLE = ParagraphStyle("title", parent=ss["Title"], fontName="DV-Bold", fontSize=19, leading=23)
SUB = ParagraphStyle("sub", parent=ss["Normal"], fontName="DV-Obl", fontSize=11,
                     textColor=colors.HexColor("#555555"), spaceAfter=2)
CAP = ParagraphStyle("cap", parent=ss["Normal"], fontName="DV-Obl", fontSize=9,
                     textColor=colors.HexColor("#555555"), leading=12, alignment=TA_LEFT, spaceBefore=3, spaceAfter=10)
BOXT = ParagraphStyle("boxt", parent=BODY, fontSize=10.5, leading=15, spaceAfter=0)

def P(t): return Paragraph(t, BODY)
def bullets(items):
    return ListFlowable([ListItem(Paragraph(x, BODY), leftIndent=10, value="•") for x in items],
                        bulletType="bullet", bulletFontName="DV", bulletFontSize=9,
                        leftIndent=14, spaceAfter=7)
def img(path, w):
    from PIL import Image as PImage
    iw, ih = PImage.open(path).size
    return Image(path, width=w, height=w*ih/iw)

S = []
S += [Paragraph("Transfer Learning from Fermi to Swift for GRB Redshift Prediction", TITLE)]
S += [Spacer(1, 4), Paragraph("A feasibility note: what we are transferring, how it is done, and what we found", SUB)]
S += [Paragraph("Prepared for the Dainotti GRB collaboration &middot; October 2026", SUB)]
S += [Spacer(1, 8), HRFlowable(width="100%", color=colors.HexColor("#cccccc"))]

# bottom line up front
box = Table([[Paragraph("<b>Bottom line.</b> We tested whether prompt gamma-ray emission information from Fermi-GBM "
    "&mdash; in particular the spectral peak energy E<sub>peak</sub> &mdash; can improve the Swift X-ray redshift "
    "estimator, by “transferring” it from the Fermi domain into the Swift domain using the 67 bursts observed "
    "by both instruments. We tried three ways of doing the transfer. <b>None improves redshift prediction.</b> The "
    "reason is physical, not a coding choice: the <i>observed</i> E<sub>peak</sub> is essentially uncorrelated with "
    "redshift (r = 0.02 across 163 bursts), so there is no redshift signal to transfer in the first place.", BOXT)]],
    colWidths=[6.7*inch])
box.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#eef3f9")),
    ("BOX",(0,0),(-1,-1),0.6,colors.HexColor("#9ab6d6")),("LEFTPADDING",(0,0),(-1,-1),10),
    ("RIGHTPADDING",(0,0),(-1,-1),10),("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8)]))
S += [Spacer(1,8), box]

S += [Paragraph("1.&nbsp;&nbsp;What are we trying to do?", H1)]
S += [P("The existing estimator predicts a gamma-ray burst&rsquo;s redshift from Swift X-ray “plateau” "
    "features &mdash; properties of the slowly-declining X-ray afterglow (the plateau flux and end-time, the decay "
    "slope, the absorbing column, and so on). It reaches a correlation of about r&nbsp;=&nbsp;0.53 between predicted "
    "and true redshift on the 281-burst sample.")]
S += [P("The proposal was to give the model <i>new</i> information from a different phase of the burst: the "
    "<b>prompt</b> gamma-ray emission, which Fermi&rsquo;s Gamma-ray Burst Monitor (GBM) records. In particular, the "
    "Fermi prompt spectrum has a well-defined <b>peak energy, E<sub>peak</sub></b> &mdash; the photon energy at which "
    "the burst radiates most of its power. E<sub>peak</sub> is physically central to GRBs (it appears in the famous "
    "Amati and Yonetoku relations), so it is a natural thing to try to fold in.")]
S += [P("The obstacle: almost none of the 281 Swift GRBs have Fermi measurements. We cannot simply paste an "
    "E<sub>peak</sub> column onto the Swift table. We need a way to carry E<sub>peak</sub> from the instrument that "
    "measures it (Fermi) into the sample that does not (Swift). <b>That bridge is what “transfer learning” "
    "refers to here.</b>")]

S += [Paragraph("2.&nbsp;&nbsp;What is &ldquo;transfer learning&rdquo; in this context? (question a)", H1)]
S += [P("In general, <b>transfer learning</b> means reusing knowledge gained in one setting &mdash; the "
    "<i>source domain</i> &mdash; to help a task in a different but related setting &mdash; the <i>target domain</i>. "
    "It is most useful when the target has limited data, or is missing a measurement that the source has. In our case:")]
S += [bullets([
    "<b>Source domain:</b> Fermi-GBM, the prompt gamma-ray phase. It measures E<sub>peak</sub> (and other prompt spectral quantities).",
    "<b>Target domain:</b> Swift, the X-ray afterglow phase. It has the 281-burst redshift sample, but no E<sub>peak</sub>.",
    "<b>The bridge:</b> 67 GRBs observed by <i>both</i> instruments. For these we have the Swift features, the Fermi E<sub>peak</sub>, and a measured redshift all at once. They are what lets us learn how the two domains relate.",
])]
S += [P("So, to answer the question directly: <b>what we are transfer-learning is E<sub>peak</sub></b> (and, more "
    "broadly, prompt-emission information). We are trying to take a quantity that lives in the Fermi world and make "
    "it usable in the Swift world, using the 67 common bursts to learn the connection. This mirrors the group&rsquo;s "
    "earlier <i>optical&nbsp;&rarr;&nbsp;X-ray</i> transfer, which brought optical-plateau information into the X-ray "
    "analysis; here the source is prompt emission instead of optical.")]

S += [Paragraph("3.&nbsp;&nbsp;The data", H1)]
S += [P("From Maria&rsquo;s spreadsheet (Sheet&nbsp;4), 67 GRBs are common to the Swift X-ray plateau sample and the "
    "Fermi-GBM catalog <i>and</i> have a measured redshift. These are the bridge bursts. For each one we retrieved the "
    "Fermi prompt spectral fit from the official Fermi-GBM Burst Catalog (the HEASARC <font face=\"DVMono\">fermigbrst</font> "
    "table), using the time-integrated Comptonized model (FLNC&nbsp;COMP), as specified: E<sub>peak</sub>, the prompt "
    "photon index, the prompt fluence and flux, each with its uncertainty.")]
S += [P("Of the 67, <b>60 have a well-constrained E<sub>peak</sub></b>; the rest have an E<sub>peak</sub> that the fit "
    "cannot pin down (their spectra are too featureless within the GBM band, so the “peak” runs off to "
    "implausible energies). We trained on the 60 good ones. For reference, a reminder of what E<sub>peak</sub> is: the "
    "energy at which the prompt spectrum peaks in nu-F-nu; the observed values here span roughly 9 to 2000 keV.")]

S += [Paragraph("4.&nbsp;&nbsp;How is the transfer done? (question b)", H1)]
S += [P("We tried three versions, from the simplest sanity check to the full “combine the samples” idea.")]
S += [Paragraph("4a.&nbsp;Direct check (on the bridge bursts only)", ParagraphStyle("h2",parent=H1,fontSize=11,textColor=colors.HexColor("#333333")))]
S += [P("The simplest test. On the 60 bridge bursts &mdash; where we have both Swift features <i>and</i> a measured "
    "E<sub>peak</sub> &mdash; we just add E<sub>peak</sub> as an extra column to the redshift model and ask whether "
    "prediction improves. This is not something we could deploy (it only works where E<sub>peak</sub> is actually "
    "measured), but it is the essential prerequisite: <i>if measured E<sub>peak</sub> does not help even here, nothing "
    "fancier will.</i>")]
S += [Paragraph("4b.&nbsp;Transfer map / imputation (the main version)", ParagraphStyle("h2b",parent=H1,fontSize=11,textColor=colors.HexColor("#333333")))]
S += [P("This is the literal “obtain E<sub>peak</sub> via transfer learning” idea. In three steps: "
    "(1)&nbsp;on the 60 bridge bursts, learn a model that predicts E<sub>peak</sub> from the Swift observables; "
    "(2)&nbsp;apply that model to all 281 Swift bursts to <i>impute</i> an E<sub>peak</sub> for each, even though they "
    "have no Fermi data; (3)&nbsp;add the imputed E<sub>peak</sub> to the redshift model and check whether prediction "
    "improves.")]
S += [Paragraph("4c.&nbsp;Pooled / domain adaptation (combine the samples)", ParagraphStyle("h2c",parent=H1,fontSize=11,textColor=colors.HexColor("#333333")))]
S += [P("Rather than imputing one variable, this combines the two <i>samples</i> into one larger training set. Because "
    "the two instruments do not measure the same things, we can only pool on the variables they <i>share</i> &mdash; the "
    "prompt photon index, the duration T90, the fluence, and the peak flux &mdash; plus a 0/1 “instrument” "
    "indicator (a <i>dummy variable</i>) so the model can absorb systematic differences between Swift and Fermi. The "
    "hope is that extra training data, from a related instrument, sharpens the predictions.")]

S += [Paragraph("5.&nbsp;&nbsp;What we found", H1)]
S += [P("All three approaches give essentially no improvement. The clearest single comparison is the real estimator "
    "(the winning GLM, with its squared and interaction terms, on the full feature set) with versus without the "
    "transfer-learned E<sub>peak</sub>:")]
S += [img("redshift_pred_with_without_epeak.png", 6.6*inch)]
S += [Paragraph("Figure&nbsp;1. Predicted vs. true redshift on the 281 Swift bursts (honest cross-validation). Left: the "
    "estimator as it stands. Right: the same estimator after adding the transfer-learned E<sub>peak</sub>. The two are "
    "point-for-point indistinguishable (r = 0.531 vs. 0.529).", CAP)]
S += [P("To avoid a confusion that is easy to fall into, here is the full ladder of numbers, with each row labelled by "
    "<i>which model</i> and <i>which features</i> it uses. The headline 0.53 is the winning GLM on the full feature "
    "set; the lower numbers further down are not failures of the transfer &mdash; they come from deliberately using a "
    "simpler model or a smaller feature set in a particular test.")]

tbl = Table([
    ["Setup (281-burst sample, honest CV)", "r"],
    ["Winning GLM (squares + interactions), full plateau features  ← the headline", "0.53"],
    ["Winning GLM + transfer-learned Epeak   (test 4b)", "0.52"],
    ["Plain linear model, full 7 plateau features", "0.46"],
    ["Plain linear model, only the 4 shared variables   (basis for test 4c)", "0.29"],
], colWidths=[5.7*inch, 0.9*inch])
tbl.setStyle(TableStyle([
    ("FONT",(0,0),(-1,0),"DV-Bold",10),("FONT",(0,1),(-1,-1),"DV",9.5),
    ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#14365c")),("TEXTCOLOR",(0,0),(-1,0),colors.white),
    ("ALIGN",(1,0),(1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, colors.HexColor("#f2f5f9")]),
    ("LINEBELOW",(0,0),(-1,-1),0.4,colors.HexColor("#cccccc")),
    ("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5),("LEFTPADDING",(0,0),(-1,-1),7)]))
S += [tbl, Spacer(1,4)]
S += [Paragraph("Table&nbsp;1. The two steps down from 0.53 to 0.29 are (i) using a plain linear model instead of the "
    "GLM (&minus;0.07, the GLM&rsquo;s nonlinear terms carry real signal) and (ii) restricting to the four variables "
    "Fermi and Swift share, which discards the strong plateau features (&minus;0.17). The pooled test (4c) necessarily "
    "lives in that bottom, shared-variable world.", CAP)]
S += [P("The precise results of each approach: <b>4a</b> (measured E<sub>peak</sub> on the 60 bridge bursts): change "
    "in r of about &minus;0.04 to 0, buried in noise. <b>4b</b> (transfer-learned E<sub>peak</sub>, winning GLM on "
    "281): 0.518&nbsp;&rarr;&nbsp;0.515, and a nonlinear random-forest check gives the same null. <b>4c</b> (pooled, "
    "shared variables + dummy): 0.287&nbsp;&rarr;&nbsp;0.275 &mdash; pooling slightly <i>hurts</i>.")]

S += [Paragraph("6.&nbsp;&nbsp;Why it does not help (the important part)", H1)]
S += [P("Three reasons, each worth understanding because they are the real conclusions:")]
S += [Paragraph("Reason 1 &mdash; there is no redshift signal in observed E<sub>peak</sub> to begin with.", ParagraphStyle("r1",parent=BODY,fontName="DV-Bold",spaceAfter=3))]
S += [P("We checked this directly on a much larger, cleaner sample: 163 Fermi GRBs with measured redshift. The "
    "observed E<sub>peak</sub> is <b>uncorrelated with redshift</b> (Pearson r = 0.02, p = 0.8).")]
S += [img("Epeak_vs_redshift_for_Maria.png", 4.5*inch)]
S += [Paragraph("Figure&nbsp;2. Observed E<sub>peak</sub> vs. redshift for 163 Fermi GRBs. The relationship is flat.", CAP)]
S += [P("The physics behind this is clean. The physically meaningful quantity is the <i>rest-frame</i> peak energy, "
    "E<sub>peak,rest</sub> = E<sub>peak,obs</sub> &times; (1+z). The Amati and Yonetoku relations link "
    "E<sub>peak,rest</sub> to a burst&rsquo;s energy or luminosity &mdash; but using them requires already knowing z. "
    "In the <i>observer</i> frame, the large intrinsic spread in E<sub>peak</sub> from burst to burst swamps the "
    "(1+z) stretching, so the observed E<sub>peak</sub> tells you almost nothing about redshift. This one fact "
    "explains every null result below it: <b>you cannot transfer information that is not there.</b>")]
S += [Paragraph("Reason 2 &mdash; the imputed E<sub>peak</sub> is redundant with the Swift features (a circularity).", ParagraphStyle("r2",parent=BODY,fontName="DV-Bold",spaceAfter=3))]
S += [P("In the transfer-map approach (4b), the imputed E<sub>peak</sub> is <i>computed from</i> the Swift features. "
    "Feeding a model a new column that is just a function of columns it already has gives it no new information. We "
    "confirmed this concretely: the imputed E<sub>peak</sub> can be reconstructed from the Swift features with r = 1.0, "
    "and adding it changes the GLM by &minus;0.004 and the random forest by &minus;0.003. Importantly, this holds for "
    "<i>any</i> model &mdash; linear, GLM, or a nonlinear ensemble like the SuperLearner &mdash; because it is about "
    "information content, not the model&rsquo;s flexibility.")]
S += [Paragraph("Reason 3 &mdash; the pooled approach can only use the weak shared variables.", ParagraphStyle("r3",parent=BODY,fontName="DV-Bold",spaceAfter=3))]
S += [P("To pool with Fermi you may only use variables both instruments measure. But Swift&rsquo;s predictive power "
    "lives in the afterglow <i>plateau</i> features (the plateau flux and time), which Fermi does not have. Dropping to "
    "the shared prompt-type variables alone costs you most of the signal (0.46&nbsp;&rarr;&nbsp;0.29) before any "
    "transfer even happens. Worse, the shared variables are not quite the same physical quantities across instruments "
    "(Swift&rsquo;s photon index is from the X-ray <i>afterglow</i>; Fermi&rsquo;s is from the <i>prompt</i> gamma "
    "rays), so their relationships to redshift differ and the pooled data pulls the fit slightly off.")]

S += [Paragraph("7.&nbsp;&nbsp;Caveats, in fairness", H1)]
S += [bullets([
    "<b>Small bridge sample.</b> The direct test (4a) uses only 60 bursts and is genuinely underpowered &mdash; which is exactly why we leaned on the 163-burst correlation (Figure&nbsp;2), which is well-powered and model-free.",
    "<b>Read the p-values loosely.</b> The significance values from repeated cross-validation are anti-conservative (the resamples reuse the same bursts), so they indicate that an effect&rsquo;s <i>direction</i> is stable, not that it is strongly significant. The effect <i>sizes</i> are tiny regardless.",
    "<b>Model choice.</b> We used the time-integrated Comptonized E<sub>peak</sub>, as specified. The power-law same-model variant is a possible follow-up, but it would not change the central fact that observed E<sub>peak</sub> is uncorrelated with z.",
    "<b>“No detectable benefit” is not “mathematically impossible.”</b> But the observer-vs-rest-frame argument makes a benefit from observed E<sub>peak</sub> specifically very unlikely.",
])]

S += [Paragraph("8.&nbsp;&nbsp;What would change the answer", H1)]
S += [P("A real gain would require a prompt-emission quantity that (a) Swift does not already have and (b) correlates "
    "with redshift <i>in the observer frame</i>. Observed E<sub>peak</sub> fails (b). The genuinely useful "
    "E<sub>peak</sub> relations are rest-frame, and using them needs z &mdash; which is circular for a redshift "
    "predictor. The other prompt variables we looked at (peak flux, fluence) were also flat-to-null, though they could "
    "be examined more carefully. If a substantially larger common sample became available, the direct test could be "
    "repeated with more statistical power &mdash; but the flat 163-burst correlation already argues strongly against a "
    "payoff.")]
S += [Spacer(1,6), HRFlowable(width="100%", color=colors.HexColor("#cccccc"))]
S += [P("<i>Data and figures underlying this note: the 67-burst bridge set (common67_transfer.csv), the 163-burst "
    "E<sub>peak</sub>&ndash;redshift set (epeak_z_163.csv), and the two figures above, are in the project&rsquo;s "
    "new_data folder.</i>")]

doc.build(S)
print("built transfer_learning_feasibility.pdf")
