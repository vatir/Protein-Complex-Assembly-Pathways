# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:48:12 2017

@author: koanb
"""
import numpy as np
import ODERun
from collections import OrderedDict

Base = ODERun.ODEFunc()
Base.GetSDataFromFiles()
Base.GetRDataFromFiles()
Base.GenCoef(1.0,1.0,1.0,1.0,BondCoef=1.0, Keff_Ident=True)

with open("Equations.tex","w") as File:
    Equation1 = ""
    Equation1 += r"\begin{equation*}"
#    Equation1 += r"k^{eff}_{i,j} = \left(\frac{k_p}{\alpha}\right)k_{d,1}^i k_{d,2}^j \mathrm{e}^{(1-i - j)9/{0.6}}"
    Equation1 += r"k^{eff}_{i,j} = \alpha \kon \kdSub{1}^i \kdSub{2}^j \mathrm{e}^{(i+j-1)\Delta G_p^0/RT}"
    Equation1 += r"\end{equation*}"+"\n"

    Equation1 += r"\begin{equation*}"
    Equation1 += r"\alpha = c_0^{-i - j}"
    Equation1 += r"\end{equation*}"+"\n"
    
    Equation1 += r"\begin{equation*}"
#    Equation1 += r"k_p = \frac{1}{2} * 10^6 \mathrm{\,M^{-1}s^{-1}}"
    Equation1 += r"k_p = \orderm{6}\ \ \konunit"
    Equation1 += r"\end{equation*}"+"\n"
    File.write(Equation1)
    File.write(r"\begin{align*}")
    
    for ID, Factors in enumerate(Base.Coef1D):
        MaxLength = 210
        ResetCount = 1

        Current = r"\DCDtNoSpacing{\species{" + "{}".format(ID+1) +r"}} &= \kon{}("
        print "Species: {}".format(ID+1)
        BiFactors = Base.Coef2D[ID]
        First = True
        for i,j in np.transpose(BiFactors.nonzero()):
            if (not First) and BiFactors[i,j] > 0:
                Current += "+"
            First = False
            if BiFactors[i,j] == 1:
                pass
            elif BiFactors[i,j] == -1:
                Current += "-"
            else:
                if BiFactors[i,j] == int(BiFactors[i,j]):
                    Current += "{}".format(int(BiFactors[i,j]))
                    print BiFactors[i,j]
                elif BiFactors[i,j] == 0.5:
                    Current += r"\frac{1}{2}"
                elif BiFactors[i,j] == 1.5:
                    Current += r"\frac{3}{2}"
                elif BiFactors[i,j] == 0.0:
                    pass
                else:
                    print "Error Bad Value Found!!!!!!!!!!!!!!!!!! : {}".format(BiFactors[i,j])
                 
            if i == j:
                Current += "\Conc{\species{"+"{}".format(i+1)+"}}^2"
            else:
                Current += "\Conc{\species{"+"{}".format(i+1)+"}}"+"\Conc{\species{"+"{}".format(j+1)+"}}"
            if len(Current) > MaxLength*ResetCount+50:
                ResetCount += 1
                Current += r"\\"+"\n"+r" & \quad "
                
        if Current[-12:] == r"\\"+"\n"+r" & \quad ":
            print "--------------"
            print Current
            Current = Current[:-12] + r") \\"+"\n"+r" & \quad "
            print Current
            print "--------------"
        else:
            Current += r")"
        MaxLength = 250
        ReacCoef = OrderedDict()
        for R in Base.Reactions:
            CurrentTerm = 0.0
            if R[0] == ID:
                CurrentTerm += R[6]
            if R[1] == ID:
                CurrentTerm += R[6]
            if R[3] == ID:
                CurrentTerm -= R[6]
            if CurrentTerm != 0.0:
                try:
                    ReacCoef[R[3]]
                except:
                    ReacCoef[R[3]] = OrderedDict()
                try:
                    ReacCoef[R[3]][(R[4],R[5])] += CurrentTerm
                except:
                    ReacCoef[R[3]][(R[4],R[5])] = CurrentTerm
            
        for Item in ReacCoef.keys():
            for BondSet in ReacCoef[Item].keys():
                print "{} : {} : {}".format(Item, BondSet, ReacCoef[Item][BondSet])
                CurrentTerm = ReacCoef[Item][BondSet]
                if CurrentTerm != 0.0:
                    if CurrentTerm > 1.0: 
                        Current += "+"
                    if CurrentTerm == 1: 
                        Current += "+"
                    elif CurrentTerm == -1.0: 
                        Current += "-"
                    elif (CurrentTerm % 1) == 0:
                        CurrentTerm = int(CurrentTerm)
                        Current += "{}".format(CurrentTerm)
                    else:
                        print CurrentTerm
                        Current += "{}".format(CurrentTerm)
                    
                    Current += "\Conc{\species{"+"{}".format(Item+1)+"}}"+r"k^{eff}_{"+"{},{}".format(BondSet[0],BondSet[1])+r"}"
                if len(Current) > MaxLength*ResetCount-25:
                    ResetCount += 1
                    Current += r"\\"+"\n"+r" & \quad "
#                Total += CurrentTerm
#            if Total != Factor:
#                print "Mismatch: {} {} {}".format(CurrentID, Total, Factor)
        if Current[-6:] == r"\quad ":
            Current = Current[:-12]
        Current += r"-\delta"+r"\Conc{\species{"+"{}".format(ID+1)+r"}}"
        if ID == 0:
            Current += "+Q"
        Current += r"\\[0.5em]"
        Current += "\n"
        print Current
        File.write(Current)
    File.write(r"\end{align*}")
    
"""
Notes:
Synth Rate for A0 equilibrium: A0*Delta
Dilution effect for all species: Delta

Keff * Combinitorial Factor * Conc1 * Conc2

"""
