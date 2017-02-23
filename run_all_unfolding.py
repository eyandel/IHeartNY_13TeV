import subprocess
import sys

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--closureTests', metavar='F', action='store_true',
                  default=False,
                  dest='closureTests',
                  help='Run closure test')

parser.add_option('--toUnfold', metavar='F', type='string', action='store',
                  default='pt',
                  dest='toUnfold',
                  help='Distribution to unfold (pt or y)')

(options, args) = parser.parse_args()
argv = []

if options.closureTests:
    path = [
        ## closure tests
        #"python unfoldTopPt.py --closureTest --toUnfold="+options.toUnfold+" --usePost",
        #"python unfoldTopPt.py --closureTest --lepType=ele --toUnfold="+options.toUnfold+" --usePost",
        #"python closureTest.py --lepType=muon --type=full --toy=nom",
        #"python closureTest.py --lepType=ele --type=full --toy=nom",
        #"python closureTest.py --lepType=muon --type=full --toy=up",
        #"python closureTest.py --lepType=ele --type=full --toy=up",
        #"python closureTest.py --lepType=muon --type=full --toy=dn",
        #"python closureTest.py --lepType=ele --type=full --toy=dn",
        "python closureTest.py --lepType=muon --type=half --toy=nom",
        "python closureTest.py --lepType=ele --type=half --toy=nom",
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --regType=LCurve",
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --regType=LCurve",
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --regType=LCurve",
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --regType=LCurve",
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --regType=LCurve",
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --regType=LCurve"
    ]
    
else :
    path = [
        "python unfoldTopPt.py --lepType=ele --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --bkgSyst=Up --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --bkgSyst=Down --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=lepUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=lepDown --toUnfold="+options.toUnfold+" --usePost",        
        "python unfoldTopPt.py --lepType=ele --systVariation=puUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=puDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=JERUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=JERDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=JECUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=JECDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=BTagUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=BTagDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=TopTagUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=TopTagDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=PDFUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=PDFDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=Q2Up --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=Q2Down --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=ASUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --lepType=ele --systVariation=ASDown --toUnfold="+options.toUnfold+" --usePost",

        "python unfoldTopPt.py --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --bkgSyst=Up --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --bkgSyst=Down --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=lepUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=lepDown --toUnfold="+options.toUnfold+" --usePost",        
        "python unfoldTopPt.py --systVariation=puUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=puDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=JERUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=JERDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=JECUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=JECDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=BTagUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=BTagDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=TopTagUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=TopTagDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=PDFUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=PDFDown --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=Q2Up --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=Q2Down --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=ASUp --toUnfold="+options.toUnfold+" --usePost",
        "python unfoldTopPt.py --systVariation=ASDown --toUnfold="+options.toUnfold+" --usePost"
    ]

## run actual unfolding
for s in path :
    print s
    subprocess.call( [s, ""], shell=True )
    
    
