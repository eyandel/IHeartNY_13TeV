import subprocess
import sys

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--closureTest', metavar='F', action='store_true',
                  default=False,
                  dest='closureTest',
                  help='Run closure test')

(options, args) = parser.parse_args()
argv = []

if options.closureTest:
    path = [
        ## closure tests
        "python unfoldTopPt.py --closureTest",
        "python unfoldTopPt.py --closureTest --lepType=ele",
    ]
else:
    path = [
        ## unfolding, combining the 1 top-tag, 1 b-tag and 1 top-tag, 0 b-tag regions
        "python unfoldTopPt.py --addNoBtag --lepType=ele",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=puUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=puDown",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=lepUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=lepDown",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=JERUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=JERDown",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=JECUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=JECDown",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=BTagUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=BTagDown",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=TopTagUp",
        "python unfoldTopPt.py --addNoBtag --lepType=ele --systVariation=TopTagDown",
        
        "python unfoldTopPt.py --addNoBtag",
        "python unfoldTopPt.py --addNoBtag --systVariation=puUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=puDown",
        "python unfoldTopPt.py --addNoBtag --systVariation=lepUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=lepDown",
        "python unfoldTopPt.py --addNoBtag --systVariation=JERUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=JERDown",
        "python unfoldTopPt.py --addNoBtag --systVariation=JECUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=JECDown",
        "python unfoldTopPt.py --addNoBtag --systVariation=BTagUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=BTagDown",
        "python unfoldTopPt.py --addNoBtag --systVariation=TopTagUp",
        "python unfoldTopPt.py --addNoBtag --systVariation=TopTagDown",
        
    ]

## run actual unfolding
for s in path :
    print s
    subprocess.call( [s, ""], shell=True )
    
    
