#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time

import ROOT

from TMVAhelper import *

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetPalette(1);
#ROOT.gROOT.ProcessLine(".L PerformanceFinal.C");

#from ROOT import myutils
ROOT.gROOT.ProcessLine(".L myutils.C++");
ROOT.gSystem.Load("myUtils_C.so");
#ROOT.gSystem.Load("PerformanceFinal_C.so");
from ROOT import myutils


############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--type',action="store",type="string",dest="type",default="gg")
parser.add_option('--ptbin',action="store",type="int",dest="ptbin",default=500)

parser.add_option('--do2D',action="store_true",dest="do2D",default=False,help='making correlations, do MVA')
parser.add_option('--train',action="store_true",dest="train",default=False,help='do training')
parser.add_option('--tmvaname',action="store",type="string",dest="tmvaname",default="qg")

(options, args) = parser.parse_args()

############################################################

# global vars
type = options.type;
bin = options.ptbin;
binp1 = options.ptbin + 100;

directory = "plots/figs_"+options.tmvaname+"_bin"+str(bin);
if not os.path.exists(directory):
    os.makedirs(directory)
############################################################

def makeCanvas(hists, names, canname, normalize=False):
    
    #bin = options.ptbin;
    #directory = "figs_bin"+str(bin);
    colors = [2,4,1,6,7];
    
    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        hists[i].SetLineColor(colors[i]);
        hists[i].SetLineWidth(2);
        if normalize: hists[i].Scale(1./hists[i].Integral());
        leg.AddEntry(hists[i], names[i], "l");

    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();

    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    can.SaveAs(directory+"/"+canname+".eps");
    can.SaveAs(directory+"/"+canname+".png");

    for hist in hists:
        hist.Scale(1./hist.Integral());

def makeCanvas2D(hist0, hist1, types):
    
    hist1.SetLineColor(2);
    
    leg = ROOT.TLegend(0.8,0.8,0.9,0.9);
    leg.SetBorderSize(1);
    leg.SetFillStyle(1);
    leg.AddEntry(hist0, types[0], "");

    leg1 = ROOT.TLegend(0.8,0.8,0.9,0.9);
    leg1.SetBorderSize(1);
    leg1.SetFillStyle(1);
    leg1.AddEntry(hist1, types[1], "");

    hist0.Scale(1./hist0.Integral());
    hist1.Scale(1./hist1.Integral());
    max = 0;
    if hist0.GetMaximum() < hist1.GetMaximum(): max = hist1.GetMaximum();
    if hist0.GetMaximum() > hist1.GetMaximum(): max = hist0.GetMaximum();

    can = ROOT.TCanvas("canOS"+hist0.GetName(),"canOS"+hist0.GetName(),1600,800);
    #hist0.SetMaximum( max*1.1 );
    #hist1.SetMaximum( max*1.1 );
    can.Divide(2,1);
    can.cd(1);
    hist0.Draw("COLZ");
    leg.Draw();
    can.cd(2);
    hist1.Draw("COLZ");
    leg1.Draw();
    can.SaveAs(directory+"/"+hist0.GetName()+"_onSame.eps");
    can.SaveAs(directory+"/"+hist0.GetName()+"_onSame.png");

def makeSingleHistCanvas( h_2D ):
    
    #directory = "figs_bin"+str(bin);

    cant = ROOT.TCanvas("cant","cant",800,800);
    h_2D.Draw("colz");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".png");
    cant.SaveAs(directory+"/"+h_2D.GetName()+".eps");

def makeROCFromHisto(hsig,hbkg,LtoR):

    nbins = hsig.GetNbinsX();
    binsize = hsig.GetBinWidth(1);
    lowedge = hsig.GetBinLowEdge(1);

    #print "lowedge: ",lowedge

    hsigIntegral = hsig.Integral();
    hbkgIntegral = hbkg.Integral();

    xval = array('d', [])
    yval = array('d', [])
    ctr = 0;
    effBkgPrev = -9999;
    for i in range(1,nbins+1):

        effBkg = 0;
        effSig = 0;

        if LtoR: effBkg = hbkg.Integral( i, nbins )/hbkgIntegral;
        else: effBkg = hbkg.Integral( 1, i )/hbkgIntegral;

        if LtoR: effSig = hsig.Integral( i, nbins )/hsigIntegral;
        else: effSig = hsig.Integral( 1, i )/hsigIntegral;

        #if not effBkg == 0.: print "cut: ",(lowedge+(i-1)*binsize),"effBkg: ", 1/effBkg, ", effSig: ", effSig;

        xval.append( effSig );
        yval.append( effBkg );

        #effBkgPrev = effBkg;
        ctr = ctr + 1;

    #print nbins, "and ", ctr
    tg = ROOT.TGraph( nbins, xval, yval );
    tg.SetName( "tg"+hsig.GetName() );
    return tg;

def SetLineColorStyleWidth(tg, color, style, width):
    tg.SetLineColor( color );
    tg.SetLineStyle( style );
    tg.SetLineWidth( width );
    
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if __name__ == '__main__':

    ROOT.TMVA.Tools.Instance();

    types = [];
    if options.tmvaname == "qg": types = ["gg","qq"];
    if options.tmvaname == "Wg": types = ["WW","gg"];
    if options.tmvaname == "Wq": types = ["WW","qq"];
    
    files = [];
    if options.tmvaname == "Wq" or options.tmvaname == "Wg": files.append("data/boost2013-jhu452-lhc8-wlwl-pt0"+str(bin)+"-0"+str(binp1)+".root");
    if options.tmvaname == "qg" or options.tmvaname == "Wg": files.append("data/boost2013-pythia81-tune4c-lhc8-gg-pt0"+str(bin)+"-0"+str(binp1)+".root");
    if options.tmvaname == "Wq" or options.tmvaname == "qg": files.append("data/boost2013-pythia81-tune4c-lhc8-qq-pt0"+str(bin)+"-0"+str(binp1)+".root");

    trainingFiles = [];
    if options.tmvaname == "Wq" or options.tmvaname == "Wg": files.append("data/boost2013-jhu452-lhc8-wlwl-pt0"+str(bin)+"-0"+str(binp1)+"_xcheck.root");
    if options.tmvaname == "qg" or options.tmvaname == "Wg": files.append("data/boost2013-pythia81-tune4c-lhc8-gg-pt0"+str(bin)+"-0"+str(binp1)+"_xcheck.root");
    if options.tmvaname == "Wq" or options.tmvaname == "qg": files.append("data/boost2013-pythia81-tune4c-lhc8-qq-pt0"+str(bin)+"-0"+str(binp1)+"_xcheck.root");

    obs_nb = 100;
    mass_obs_nb = 1000;

    #######-----------------------------------------------
    ####### setup
    #######-----------------------------------------------

    ## =====================================================================================
    ## 1D setup
    ## =====================================================================================
    
    h_pt1 = [];
    h_pt2 = [];
    h_eta1 = [];
    h_eta2 = [];
    h_mass1 = [];
    h_mass2 = [];
    h_tau1_b1 = [];
    h_tau1_b2 = [];
    h_tau21_b1 = [];
    h_tau21_b2 = [];
    h_c1_b0 = [];
    h_c1_b1 = [];
    h_c1_b2 = [];
    h_c2_b1 = [];
    h_c2_b2 = [];
    h_qjetVol = [];
    h_mass_trim = [];
    h_mass_mmdt = [];
    h_mass_prun = [];
    h_mass_sdb2 = [];
    h_mass_sdm1 = [];
    h_multiplicity = [];
    
    h_bdtall = [];
    
    for i in range(len(types)):
        
        h_pt1.append( ROOT.TH1F("h_jpt1_"+types[i],";jet 1 pT (GeV); N",50,0,2000) );
        h_pt2.append( ROOT.TH1F("h_jpt2_"+types[i],";jet 2 pT (GeV); N",50,0,2000) );
        h_eta1.append( ROOT.TH1F("h_jeta1_"+types[i],";jet 1 eta; N",20,-5,5) );
        h_eta2.append( ROOT.TH1F("h_jeta2_"+types[i],";jet 2 eta; N",20,-5,5) );
        h_mass1.append( ROOT.TH1F("h_jmass1_"+types[i],";jet 1 mass; N",mass_obs_nb,0,200) );
        h_mass2.append( ROOT.TH1F("h_jmass2_"+types[i],";jet 2 mass; N",obs_nb,0,200) );
        
        h_tau1_b1.append( ROOT.TH1F("h_jtau1_b1"+types[i],";tau1 (beta=1); N",obs_nb,0,0.5) );
        h_tau1_b2.append( ROOT.TH1F("h_jtau1_b2"+types[i],";tau1 (beta=2); N",obs_nb,0,0.5) );
        h_tau21_b1.append( ROOT.TH1F("h_jtau21_b1"+types[i],";tau21 (beta=1); N",obs_nb,0,1.0) );
        h_tau21_b2.append( ROOT.TH1F("h_jtau21_b2"+types[i],";tau21 (beta=2); N",obs_nb,0,1.0) );
        h_c1_b0.append( ROOT.TH1F("h_jc1_b0"+types[i],";c1 (beta=0); N",obs_nb,0.2,0.5) );
        h_c1_b1.append( ROOT.TH1F("h_jc1_b1"+types[i],";c1 (beta=1); N",obs_nb,0,0.3) );
        h_c1_b2.append( ROOT.TH1F("h_jc1_b2"+types[i],";c1 (beta=2); N",obs_nb,0,0.3) );
        h_c2_b1.append( ROOT.TH1F("h_jc2_b1"+types[i],";c2 (beta=1); N",obs_nb,0,1.0) );
        h_c2_b2.append( ROOT.TH1F("h_jc2_b2"+types[i],";c2 (beta=2); N",obs_nb,0,1.0) );
        
        h_qjetVol.append( ROOT.TH1F("h_j_qjetVol"+types[i],";#Gamma_{Qjet}; N",obs_nb,0,1) );
        h_mass_trim.append( ROOT.TH1F("h_j_mass_trim"+types[i],";m_{trim}; N",mass_obs_nb,0,240) );
        h_mass_mmdt.append( ROOT.TH1F("h_j_mass_mmdt"+types[i],";m_{mmdt}; N",mass_obs_nb,0,240) );
        h_mass_prun.append( ROOT.TH1F("h_j_mass_prun"+types[i],";m_{prun}; N",mass_obs_nb,0,240) );
        h_mass_sdb2.append( ROOT.TH1F("h_j_mass_sdb2"+types[i],";m_{sd}^{#beta=2}; N",mass_obs_nb,0,240) );
        h_mass_sdm1.append( ROOT.TH1F("h_mass_sdm1"+types[i],";m_{sd}^{#beta=-1}; N",mass_obs_nb,0,240) );
        h_multiplicity.append( ROOT.TH1F("h_multiplicity"+types[i],";n_{constits}; N",150,0,150) );

        h_bdtall.append( ROOT.TH1F("h_bdtall"+types[i],";bdt all; N",100,-1,1) );

    ## =====================================================================================
    ## 2D setup
    ## =====================================================================================

    if options.do2D:
        
        obs2Corr = [];
        obs2Cor_name = [];
        obs2Corr_nb = [];
        obs2Corr_lo = [];
        obs2Corr_hi = [];
        
        if options.tmvaname == "qg":
            #observables to correlate -- quark vs. gluon
            obs2Corr = ["jmass","jtau1_b1","jtau1_b2","jc1_b0","jc1_b1","jc1_b2","j_multiplicity","j_qjetVol"];
            obs2Corr_name = ["m","#tau_{1}^{#beta=1}","#tau_{1}^{#beta=2}",
                             "C_{1}^{#beta=0}","C_{1}^{#beta=1}","C_{1}^{#beta=2}",
                             "n_{constits}","#Gamma_{Qjet}"]
            #obs2Corr_nb = [40,40,40,40,40,40,40,40,40,40];
            obs2Corr_nb = [40]*8;
            obs2Corr_lo = [  0,  0,  0,0.2,  0,  0,  0,  0];
            obs2Corr_hi = [200,0.5,0.5,0.5,0.3,0.3,150,  1];

        if options.tmvaname == "Wq" or options.tmvaname == "Wg":
            #observables to correlate -- W vs. quark/gluon
            obs2Corr = ["jmass","jtau21_b1","jtau21_b2",
                        "jc2_b1","jc2_b2",
                        "j_qjetVol",
                        "j_mass_trim","j_mass_mmdt","j_mass_prun","j_mass_sdb2","j_mass_sdm1"];
            obs2Corr_name = ["m","#tau_{21}^{#beta=1}","#tau_{21}^{#beta=2}",
                             "C_{2}^{#beta=1}","C_{2}^{#beta=2}",
                             "#Gamma_{Qjet}",
                             "m_{trim}","m_{mmdt}","m_{prun}","m_{sd}^{#beta=2}","m_{sd}^{#beta=-1}"]
            obs2Corr_nb = [40]*11;
            obs2Corr_lo = [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0];
            obs2Corr_hi = [200,1.0,1.0,0.6,0.6,  1,200,200,200,200,200];

        h_bdts = [];
        h_bdts_vars = [];
        h_2D = [];
        for a in range(len(types)):
            h_2D_type = []
            h_bdts_type = [];
            for i in range(len(obs2Corr)):
                for j in range(i+1,len(obs2Corr)):
                    hname = "h2d_"+obs2Corr[i]+"_"+obs2Corr[j]+"_"+types[a];
                    htitle = ";"+obs2Corr_name[i]+";"+obs2Corr_name[j];
                    print "hname = ",hname
                    h_2D_type.append( ROOT.TH2F(hname,htitle,obs2Corr_nb[i],obs2Corr_lo[i],obs2Corr_hi[i],obs2Corr_nb[j],obs2Corr_lo[j],obs2Corr_hi[j]) );
                    
                    # bdt stuff
                    bdttitle = ";"+obs2Corr_name[i]+"+"+obs2Corr_name[j]+"; N";
                    h_bdts_type.append( ROOT.TH1F("h_bdt_"+hname,bdttitle,obs_nb,-1,1) )
                    if a == 0: h_bdts_vars.append( [obs2Corr[i],obs2Corr[j]] );
                        
            h_2D.append( h_2D_type );
            h_bdts.append( h_bdts_type );
                    
        #######-----------------------------------------------
        #######-----------------------------------------------
        ## TMVA stuff
        ## classifier
        ## inputs, list of variables, list of ntuples
        bdtClasses = [];
        bdtctr = 0;
        for i in range(len(obs2Corr)):
            for j in range(i+1,len(obs2Corr)):
                
                vars_tmp = [obs2Corr[i],obs2Corr[j]];
                tmva_tmp = TMVAhelper("tmva_"+options.tmvaname+"_"+obs2Corr[i]+"_"+obs2Corr[j],vars_tmp,files);
                bdtClasses.append(tmva_tmp);
                if options.train: bdtClasses[bdtctr].train();
                bdtClasses[bdtctr].read();            
                bdtctr += 1;

        tmva_all = TMVAhelper("tmva_"+options.tmvaname+"_all",obs2Corr,files);
        if options.train: tmva_all.train();
        tmva_all.read();
            
    ## =====================================================================================
    ## =====================================================================================
    
    #######-----------------------------------------------
    ####### loop
    #######-----------------------------------------------
    
    for a in range(len(types)):
    
        f = ROOT.TFile(files[a]);
        t = f.Get("t");
    
        for i in range(t.GetEntries()):

            if i % 1000 == 0: print "file #",a,", event #",i
            if i > 100000: break;
            
            t.GetEntry(i);
            h_pt1[a].Fill( t.jpt[0] );
            h_pt2[a].Fill( t.jpt[1] );
            h_eta1[a].Fill( t.jeta[0] );
            h_eta2[a].Fill( t.jeta[1] );
            h_mass1[a].Fill( t.jmass[0] );
            h_mass2[a].Fill( t.jmass[1] );
            
            h_tau1_b1[a].Fill( t.jtau1_b1[0] );
            h_tau1_b2[a].Fill( t.jtau1_b2[0] );            
            h_tau21_b1[a].Fill( t.jtau21_b1[0] );
            h_tau21_b2[a].Fill( t.jtau21_b2[0] );
            
            h_c1_b0[a].Fill( t.jc1_b0[0] );                        
            h_c1_b1[a].Fill( t.jc1_b1[0] );                        
            h_c1_b2[a].Fill( t.jc1_b2[0] ); 
            h_c2_b1[a].Fill( t.jc2_b1[0] );
            h_c2_b2[a].Fill( t.jc2_b2[0] );
            
            h_qjetVol[a].Fill( t.j_qjetVol[0] );
            h_mass_trim[a].Fill( t.j_mass_trim[0] );
            h_mass_mmdt[a].Fill( t.j_mass_mmdt[0] );
            h_mass_prun[a].Fill( t.j_mass_prun[0] );
            h_mass_sdb2[a].Fill( t.j_mass_sdb2[0] );
            h_mass_sdm1[a].Fill( t.j_mass_sdm1[0] );
            h_multiplicity[a].Fill( t.j_multiplicity[0] );
            
            ## =====================================================================================
            ## 2D stuff
            if options.do2D:
                #2D plots
                obsctr = 0;
                for j in range(len(obs2Corr)):
                    for k in range(j+1,len(obs2Corr)):
                        h_2D[a][obsctr].Fill( getattr(t,obs2Corr[j])[0], getattr(t,obs2Corr[k])[0] );

                        tmpList = [];
                        for l in range(len(h_bdts_vars[obsctr])): 
                            tmpList.append( getattr(t,h_bdts_vars[obsctr][l])[0] );
                        bdtOutput = bdtClasses[obsctr].evaluate( tmpList );
                        h_bdts[a][obsctr].Fill(bdtOutput);
                        obsctr += 1;
            
                tmpList = [];
                for j in range(len(obs2Corr)):
                    tmpList.append( getattr(t,obs2Corr[j])[0] );
                h_bdtall[a].Fill( tmva_all.evaluate( tmpList ) );
            ## =====================================================================================

        del f;
        del t;
                
    #######-----------------------------------------------
    ####### plot
    #######-----------------------------------------------

    ## =====================================================================================
    ## 1D stuff

    makeCanvas(h_pt1, types, "jpt1", True);
    makeCanvas(h_pt2, types, "jpt2", True);
    makeCanvas(h_eta1, types, "jeta1", True);
    makeCanvas(h_eta2, types, "jeta2", True);
    makeCanvas(h_mass1, types, "jmass1", True);
    makeCanvas(h_mass2, types, "jmass2", True);
                
    makeCanvas(h_tau1_b1, types, "h_tau1_b1", True);
    makeCanvas(h_tau1_b2, types, "h_tau1_b2", True);
    makeCanvas(h_tau21_b1, types, "h_tau21_b1", True);
    makeCanvas(h_tau21_b2, types, "h_tau21_b2", True);
    makeCanvas(h_c1_b0, types, "h_c1_b0", True);
    makeCanvas(h_c1_b1, types, "h_c1_b1", True);
    makeCanvas(h_c1_b2, types, "h_c1_b2", True);
    makeCanvas(h_c2_b1, types, "h_c2_b1", True);
    makeCanvas(h_c2_b2, types, "h_c2_b2", True);

    makeCanvas(h_qjetVol, types, "h_qjetVol", True);
    makeCanvas(h_mass_trim, types, "h_mass_trim", True);
    makeCanvas(h_mass_mmdt, types, "h_mass_mmdt", True);
    makeCanvas(h_mass_prun, types, "h_mass_prun", True);
    makeCanvas(h_mass_sdb2, types, "h_mass_sdb2", True);
    makeCanvas(h_mass_sdm1, types, "h_mass_sdm1", True);
    makeCanvas(h_multiplicity, types, "h_multiplicity", True);

    makeCanvas(h_bdtall, types, "h_bdtall", True);

    ## =====================================================================================
    ## 2D stuff
    if options.do2D:
        
        #plot bdts
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                tmpList = [];
                for a in range(len(types)): tmpList.append( h_bdts[a][obsctr] );
                makeCanvas( tmpList, types, h_bdts[0][obsctr].GetName() );
                obsctr += 1;

        for a in range(len(types)):
            #2D plots
            obsctr = 0;
            for j in range(len(obs2Corr)):
                for k in range(j+1,len(obs2Corr)):
                    makeSingleHistCanvas(h_2D[a][obsctr]);        
                    obsctr += 1;
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                makeCanvas2D(h_2D[0][obsctr],h_2D[1][obsctr],types);
                obsctr += 1;

        # for Ben's utility
        theutils = myutils()
        ngp = 25;

        # single ROCs
        roc_mass    = theutils.ROC( h_mass1[0], h_mass1[1], ngp );
        SetLineColorStyleWidth( roc_mass, 1, 1, 6 );

        roc_tau1_b1 = makeROCFromHisto( h_tau1_b1[0], h_tau1_b1[1], 1 );
        SetLineColorStyleWidth( roc_tau1_b1, 2, 9, 6 );
        roc_tau1_b2 = makeROCFromHisto( h_tau1_b2[0], h_tau1_b2[1], 1 );
        SetLineColorStyleWidth( roc_tau1_b2, 4, 9, 6 );
        roc_tau21_b1 = makeROCFromHisto( h_tau21_b1[0], h_tau21_b1[1], 0 );
        SetLineColorStyleWidth( roc_tau21_b1, 2, 9, 6 );
        roc_tau21_b2 = makeROCFromHisto( h_tau21_b2[0], h_tau21_b2[1], 0  );
        SetLineColorStyleWidth( roc_tau21_b2, 4, 9, 6 );
        
        roc_c1_b0   = makeROCFromHisto( h_c1_b0[0], h_c1_b0[1], 1 );
        SetLineColorStyleWidth( roc_c1_b0, 6, 9, 6 );
        roc_c1_b1   = makeROCFromHisto( h_c1_b1[0], h_c1_b1[1], 1 );
        SetLineColorStyleWidth( roc_c1_b1, 7, 9, 6 );
        roc_c1_b2   = makeROCFromHisto( h_c1_b2[0], h_c1_b2[1], 1 );
        SetLineColorStyleWidth( roc_c1_b2, 3, 9, 6 );
        roc_c2_b1   = makeROCFromHisto( h_c2_b1[0], h_c2_b1[1], 0 );
        SetLineColorStyleWidth( roc_c2_b1, 7, 9, 6 );
        roc_c2_b2   = makeROCFromHisto( h_c2_b2[0], h_c2_b2[1], 0 );
        SetLineColorStyleWidth( roc_c2_b2, 3, 9, 6 );

        roc_qjetVol = makeROCFromHisto( h_qjetVol[0], h_qjetVol[1], 0 );
        SetLineColorStyleWidth( roc_qjetVol, 3, 1, 6 );

        roc_mass_trim = theutils.ROC( h_mass_trim[0], h_mass_trim[1], ngp );
        SetLineColorStyleWidth( roc_mass_trim, 2, 1, 6 );
        roc_mass_mmdt = theutils.ROC( h_mass_mmdt[0], h_mass_mmdt[1], ngp );
        SetLineColorStyleWidth( roc_mass_mmdt, 15, 1, 6 );
        roc_mass_prun = theutils.ROC( h_mass_prun[0], h_mass_prun[1], ngp );
        SetLineColorStyleWidth( roc_mass_mmdt, 6, 1, 6 );
        roc_mass_sdb2 = theutils.ROC( h_mass_sdb2[0], h_mass_sdb2[1], ngp );
        SetLineColorStyleWidth( roc_mass_sdb2, 7, 1, 6 );
        roc_mass_sdm1 = theutils.ROC( h_mass_sdm1[0], h_mass_sdm1[1], ngp );
        SetLineColorStyleWidth( roc_mass_sdm1, 4, 1, 6 );

        roc_multiplicity = makeROCFromHisto( h_multiplicity[0], h_multiplicity[1], 1 );
        SetLineColorStyleWidth( roc_mass_sdm1, 1, 9, 6 );

        roc_bdtall = makeROCFromHisto( h_bdtall[0], h_bdtall[1], 1 );
        SetLineColorStyleWidth( roc_bdtall, 28, 1, 6 );


        #bdt styles...
        bdtcolors = [1,2,3,4,6,7]
        bdtstyles = [1,3,5,7,9]
        bdtwidths = [1,3];

        bdtlines = [];
        for j in range(len(bdtstyles)):
            for k in range(len(bdtcolors)):
                for m in range(len(bdtwidths)):
                    bdtlines.append( [bdtcolors[k], bdtstyles[j], bdtwidths[m]] );
        roc_bdts = [];
        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                roc_bdts.append( makeROCFromHisto( h_bdts[0][obsctr], h_bdts[1][obsctr], 1 ) );
                roc_bdts[obsctr].SetLineColor( bdtlines[obsctr][0] );
                roc_bdts[obsctr].SetLineStyle( bdtlines[obsctr][1] );
                roc_bdts[obsctr].SetLineWidth( bdtlines[obsctr][2] );
                obsctr += 1;

        ## combine all ROCs into one array...
        allROCs = [];
        allROCsNames = [];
        if options.tmvaname == "Wq" or options.tmvaname == "Wg" or options.tmvaname == "qg":
            allROCs.append( roc_mass );
            allROCsNames.append( "m" );
        if options.tmvaname == "qg":
            allROCs.append( roc_tau1_b1 );
            allROCsNames.append( "#tau_{1}^{#beta=1}" );
            allROCs.append( roc_tau1_b2 );  
            allROCsNames.append( "#tau_{1}^{#beta=2}" );
        if options.tmvaname == "Wq" or options.tmvaname == "Wg":
            allROCs.append( roc_tau21_b1 );
            allROCsNames.append( "#tau_{21}^{#beta=1}" );
            allROCs.append( roc_tau21_b2 );
            allROCsNames.append( "#tau_{21}^{#beta=2}" );
        if options.tmvaname == "qg":
            allROCs.append( roc_c1_b0 );
            allROCsNames.append( "C_{1}^{#beta=0}" );
            allROCs.append( roc_c1_b1 );  
            allROCsNames.append( "C_{1}^{#beta=1}" );
            allROCs.append( roc_c1_b2 );  
            allROCsNames.append( "C_{1}^{#beta=2}" );
        if options.tmvaname == "Wq" or options.tmvaname == "Wg":
            allROCs.append( roc_c2_b1 );
            allROCsNames.append( "C_{2}^{#beta=1}" );
            allROCs.append( roc_c2_b2 );
            allROCsNames.append( "C_{2}^{#beta=2}" );
        if options.tmvaname == "Wq" or options.tmvaname == "Wg" or options.tmvaname == "qg":
            allROCs.append( roc_qjetVol );
            allROCsNames.append( "#Gamma_{Qjet}" );
        if options.tmvaname == "Wq" or options.tmvaname == "Wg":
            allROCs.append( roc_mass_trim );
            allROCsNames.append( "m_{trim}" );
            allROCs.append( roc_mass_mmdt );
            allROCsNames.append( "m_{mmdt}" );
            allROCs.append( roc_mass_prun );  
            allROCsNames.append( "m_{prun}" );
            allROCs.append( roc_mass_sdb2 );  
            allROCsNames.append( "m_{sd}^{#beta=2}" );
            allROCs.append( roc_mass_sdm1 );
            allROCsNames.append( "m_{sd}^{#beta=-1}" );
        if options.tmvaname == "qg":
            allROCs.append( roc_multiplicity );
            allROCsNames.append( "n_{constits}" );

        obsctr = 0;
        for j in range(len(obs2Corr)):
            for k in range(j+1,len(obs2Corr)):
                allROCs.append( roc_bdts[obsctr] );  
                allROCsNames.append( obs2Corr_name[j]+"+"+obs2Corr_name[k] );
                obsctr += 1;

        allROCs.append( roc_bdtall );
        allROCsNames.append( "allvars" );

        ## --------------------
        ## start drawing
        leg = ROOT.TLegend( 0.0, 0.0, 1.0, 1.0 );
        leg.SetBorderSize( 0 );
        leg.SetFillStyle( 0 );
        #leg.SetTextSize( 0.02 );
        leg.SetNColumns(4);
        for i in range(len(allROCs)):
            leg.AddEntry( allROCs[i], allROCsNames[i], "l" );

        canRoc = ROOT.TCanvas("canRoc","canRoc",2500,1000);    
        canRoc.Divide(2,1);
        canRoc.cd(2);
        leg.Draw();
        canRoc.cd(1);
        hrl = canRoc.DrawFrame(0,1e-3,1.0,1);
        hrl.GetXaxis().SetTitle("#epsilon_{sig}");
        hrl.GetYaxis().SetTitle("#epsilon_{bkg}");
        for i in range(len(allROCs)):
            allROCs[i].Draw();
        canRoc.SaveAs(directory+"/Rocs_1D_lin.eps");
        canRoc.SaveAs(directory+"/Rocs_1D_lin.png");
        canRoc.SaveAs(directory+"/Rocs_1D_lin.pdf");
        ROOT.gPad.SetLogy();
        canRoc.SaveAs(directory+"/Rocs_1D.eps");
        canRoc.SaveAs(directory+"/Rocs_1D.png");
        canRoc.SaveAs(directory+"/Rocs_1D.pdf");

        ## ------------
        
        for a in range(len(obs2Corr)):
            
            legSingle = ROOT.TLegend( 0.0, 0.0, 1.0, 1.0 );
            legSingle.SetBorderSize( 0 );
            legSingle.SetFillStyle( 0 );
            #leg.SetTextSize( 0.06 );
            legSingle.SetNColumns(4);
            for i in range(len(allROCs)):
                if obs2Corr[a] in allROCs[i].GetName(): legSingle.AddEntry( allROCs[i], allROCsNames[i], "l" );

            canRocSingle = ROOT.TCanvas("canRocSingle","canRocSingle",2500,1000);    
            canRocSingle.Divide(2,1);
            canRocSingle.cd(1);
            hrlSingle = canRocSingle.DrawFrame(0,1e-3,1.0,1);
            hrlSingle.GetXaxis().SetTitle("#epsilon_{sig}");
            hrlSingle.GetYaxis().SetTitle("#epsilon_{bkg}");
            for i in range(len(allROCs)):
                if obs2Corr[a] in allROCs[i].GetName(): allROCs[i].Draw();
            ROOT.gPad.SetLogy();

            canRocSingle.cd(2);
            legSingle.Draw();
            canRocSingle.SaveAs(directory+"/Rocs_1D_"+obs2Corr[a]+".eps");
            canRocSingle.SaveAs(directory+"/Rocs_1D_"+obs2Corr[a]+".png");

            del legSingle
            del hrlSingle
            del canRocSingle
            
        ## ------------    
        bdtallIndex = len(allROCs)-1;
        leg1 = ROOT.TLegend( 0.0, 0.0, 1.0, 1.0 );
        leg1.SetBorderSize( 0 );
        leg1.SetFillStyle( 0 );
        #leg1.SetTextSize( 0.01 );
        leg1.SetNColumns(4);
        for i in range(len(obs2Corr)):
            leg1.AddEntry( allROCs[i], allROCsNames[i], "l" );
        leg1.AddEntry( allROCs[bdtallIndex], allROCsNames[bdtallIndex], "l" );

        canRoc1 = ROOT.TCanvas("canRoc1","canRoc1",2500,1000);    
        canRoc1.Divide(2,1);
        canRoc1.cd(2);
        leg1.Draw();
        canRoc1.cd(1);
        hrl1 = canRoc1.DrawFrame(0,1e-3,1.0,1);
        hrl1.GetXaxis().SetTitle("#epsilon_{sig}");
        hrl1.GetYaxis().SetTitle("#epsilon_{bkg}");
        for i in range(len(obs2Corr)):
            allROCs[i].Draw("l");
        allROCs[bdtallIndex].Draw("l");
        canRoc1.SaveAs(directory+"/Rocs_1D_single_lin.eps");
        canRoc1.SaveAs(directory+"/Rocs_1D_single_lin.png");
        canRoc1.SaveAs(directory+"/Rocs_1D_single_lin.pdf");
        ROOT.gPad.SetLogy();
        canRoc1.SaveAs(directory+"/Rocs_1D_single.eps");
        canRoc1.SaveAs(directory+"/Rocs_1D_single.png");        
        canRoc1.SaveAs(directory+"/Rocs_1D_single.pdf");
